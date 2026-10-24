function result = qrrt_jeong2019(s,overrides)
if nargin<2, overrides=struct; end
[o,T,result,state,t0]=planners.fq_support('setup',s,overrides,'Q');
trace=struct([]); draws=0; accepted=0; rewires=0;
reason='iteration_limit'; stop_first=strcmp(o.stop_mode,'first_solution');
if result.success && stop_first, reason='first_solution'; end
for iter=1:o.max_iter
    if result.success && stop_first, break; end
    if toc(t0)>=o.time_limit, reason='time_limit'; break; end
    result.iterations=iter;
    [sample,ok,sample_reason,count]=planners.fq_support('sample',s,o,iter,t0);
    draws=draws+count;
    if ~ok
        if strcmp(sample_reason,'sample_in_obstacle'), continue; end
        reason=sample_reason; break;
    end
    [nearest,dist]=rr.nearest(T,sample); origin=T.X(nearest,:);
    q=sample;
    if dist>o.step, q=origin+o.step*(sample-origin)/dist; end
    tr=struct('iteration',iter,'sample',sample,'nearest_id',nearest, ...
        'steer_origin',origin,'new_point',q,'accepted',false,'new_id',0, ...
        'chosen_parent_id',0,'near_ids',zeros(1,0),'ancestor_ids',zeros(1,0), ...
        'parent_candidates',zeros(1,0),'rewire_sources',zeros(1,0), ...
        'rewired',zeros(0,3));
    if o.record_trace && isempty(trace), trace=tr([]); end
    if dist==0 || rr.collision(origin,q,s)
        if o.record_trace, trace(end+1)=tr; end %#ok<AGROW>
        continue;
    end
    near=find(vecnorm(T.X-q,2,2)<=o.near_radius)';
    ancestors=ancestry(T,near,o.depth);
    candidates=unique([near ancestors],'stable');
    parent=nearest; best=T.c(nearest)+norm(q-origin);
    for v=candidates
        cost=T.c(v)+norm(T.X(v,:)-q);
        if cost+1e-12<best && ~rr.collision(T.X(v,:),q,s)
            parent=v; best=cost;
        end
    end
    [T,id]=rr.insert(T,q,parent); accepted=accepted+1;
    tr.accepted=true; tr.new_id=id; tr.chosen_parent_id=parent;
    tr.near_ids=near; tr.ancestor_ids=ancestors; tr.parent_candidates=candidates;
    [result,state]=planners.fq_support('observe',T,id,s,o,result,state,t0,iter);
    if result.success && stop_first
        reason='first_solution';
        if o.record_trace, trace(end+1)=tr; end %#ok<AGROW>
        break;
    end
    sources=[id T.p(id)]; tr.rewire_sources=sources;
    for v=near
        for source=sources
            if v==1 || v==source, continue; end
            cost=T.c(source)+norm(T.X(source,:)-T.X(v,:));
            if cost+1e-12<T.c(v) && ~rr.collision(T.X(source,:),T.X(v,:),s)
                old=T.p(v); T=rr.reparent(T,v,source);
                tr.rewired(end+1,:)=[v old source]; rewires=rewires+1;
            end
        end
    end
    [result,state]=planners.fq_support('observe',T,[],s,o,result,state,t0,iter);
    if o.record_trace, trace(end+1)=tr; end %#ok<AGROW>
end
details=struct('source_doi','10.1016/j.eswa.2019.01.032', ...
    'source_algorithms','Article p84 Algorithm 2; p85 Algorithms 4-5 and section 3.1.2', ...
    'sampling_draws',draws,'accepted_samples',accepted,'rewirings',rewires, ...
    'native_termination','Full iteration budget (Algorithm 4); wall-clock cap is an experiment safeguard.', ...
    'rewire_rule','New node and immediate parent; independent of Depth.');
result=planners.fq_support('finish',T,s,o,result,state,t0,trace,details,reason);
end

function ids=ancestry(T,vertices,depth)
ids=zeros(1,0);
for v=vertices
    for k=1:depth
        v=T.p(v);
        if v==0, break; end
        ids(end+1)=v; %#ok<AGROW>
    end
end
ids=unique(ids,'stable');
end
