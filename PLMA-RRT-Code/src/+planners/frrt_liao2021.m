function result = frrt_liao2021(s,overrides)
if nargin<2, overrides=struct; end
[o,T,result,state,t0]=planners.fq_support('setup',s,overrides,'F');
trace=struct([]); draws=0; accepted=0; created=0; rewires=0;
reason='iteration_limit'; stop_first=~strcmp(o.stop_mode,'budget');
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
        'chosen_parent_id',0,'near_ids',zeros(1,0),'reachest_id',0, ...
        'tested_ancestors',zeros(1,0),'created_id',0,'create_point',[], ...
        'dichotomy_steps',[0 0],'rewired',zeros(0,3));
    if o.record_trace && isempty(trace), trace=tr([]); end
    if dist==0 || rr.collision(origin,q,s)
        if o.record_trace, trace(end+1)=tr; end %#ok<AGROW>
        continue;
    end
    near=find(vecnorm(T.X-q,2,2)<=o.near_radius)';
    [reachest,tested]=find_reachest(T,nearest,q,s);
    [qc,steps]=create_node(T,reachest,q,s,o.dichotomy_distance);
    ids=zeros(1,0); parent=reachest;
    if ~isempty(qc)
        [T,parent]=rr.insert(T,qc,T.p(reachest));
        ids(end+1)=parent; created=created+1;
        tr.created_id=parent; tr.create_point=qc;
    end
    [T,id]=rr.insert(T,q,parent); ids(end+1)=id; accepted=accepted+1;
    tr.accepted=true; tr.new_id=id; tr.chosen_parent_id=parent;
    tr.near_ids=near; tr.reachest_id=reachest;
    tr.tested_ancestors=tested; tr.dichotomy_steps=steps;
    [result,state]=planners.fq_support('observe',T,ids,s,o,result,state,t0,iter);
    if result.success && stop_first
        reason='first_solution';
        if o.record_trace, trace(end+1)=tr; end %#ok<AGROW>
        break;
    end
    for v=near
        cost=T.c(id)+norm(T.X(id,:)-T.X(v,:));
        if v~=1 && cost+1e-12<T.c(v) && ~rr.collision(T.X(id,:),T.X(v,:),s)
            old=T.p(v); T=rr.reparent(T,v,id);
            tr.rewired(end+1,:)=[v old id]; rewires=rewires+1;
        end
    end
    [result,state]=planners.fq_support('observe',T,[],s,o,result,state,t0,iter);
    if o.record_trace, trace(end+1)=tr; end %#ok<AGROW>
end
details=struct('source_doi','10.1016/j.eswa.2021.115457', ...
    'source_algorithms','PDF p6 Algorithm 7; p7 Algorithm 8; p8 Algorithm 9; p3 Algorithm 3', ...
    'sampling_draws',draws,'accepted_samples',accepted,'created_nodes',created, ...
    'rewirings',rewires,'native_termination','First feasible solution (Algorithm 7 lines 16-18).', ...
    'budget_mode_extension',strcmp(o.stop_mode,'budget'));
result=planners.fq_support('finish',T,s,o,result,state,t0,trace,details,reason);
end

function [reachest,tested]=find_reachest(T,nearest,q,s)
reachest=nearest; tested=zeros(1,0);
while T.p(reachest)~=0
    p=T.p(reachest); tested(end+1)=p; %#ok<AGROW>
    if rr.collision(q,T.X(p,:),s), return; end
    reachest=p;
end
end

function [qc,steps]=create_node(T,reachest,q,s,tol)
qc=[]; steps=[0 0]; p=T.p(reachest);
if p==0, return; end
allow=T.X(reachest,:); forbid=T.X(p,:);
while norm(allow-forbid)>tol
    mid=(allow+forbid)/2;
    if isequal(mid,allow) || isequal(mid,forbid), break; end
    steps(1)=steps(1)+1;
    if ~rr.collision(q,mid,s), allow=mid; else, forbid=mid; end
end
forbid=q;
while norm(allow-forbid)>tol
    mid=(allow+forbid)/2;
    if isequal(mid,allow) || isequal(mid,forbid), break; end
    steps(2)=steps(2)+1;
    if ~rr.collision(mid,T.X(p,:),s), allow=mid; else, forbid=mid; end
end
if ~isequal(allow,T.X(reachest,:)) && ...
        ~rr.collision(T.X(p,:),allow,s) && ~rr.collision(allow,q,s)
    qc=allow;
end
end
