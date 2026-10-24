function result = bi_jordan2013(s,overrides)
if nargin<2, overrides=struct; end
o=options(s,overrides); D=numel(s.start); t0=tic;
if rr.collision(s.start,s.start,s) || rr.collision(s.goal,s.goal,s)
    error('jordan:Endpoints','Start and goal must be in free space.');
end
A=rr.new_tree(s.start); B=rr.new_tree(s.goal); A.active=true; B.active=true;
links=zeros(0,2); best=Inf; best_link=0;
result=struct('path',zeros(0,D),'raw_path',zeros(0,D),'treeA',A,'treeB',B, ...
    'iterations',0,'success',false,'first_time',NaN,'elapsed',0, ...
    'history',zeros(0,3),'details',struct);
first_path=zeros(0,D); first_iteration=NaN; trace=struct([]);
draws=0; accepted=0; rewirings=0; connections=0; contractions=0; pruned=0;
reason='iteration_limit';
if isequal(s.start(:)',s.goal(:)')
    links=[1 1]; [result,best,best_link,first_path,first_iteration]= ...
        observe(A,B,links,result,best,first_path,first_iteration,t0,0);
end
for iter=1:o.max_iter
    if result.success && strcmp(o.stop_mode,'first_solution'), reason='first_solution'; break; end
    if toc(t0)>=o.time_limit, reason='time_limit'; break; end
    result.iterations=iter; growA=mod(iter,2)==1;
    if growA, T=A; U=B; target=s.goal; else, T=B; U=A; target=s.start; end
    [sample,ok,sample_reason,n]=draw_sample(s,o,iter,target,t0); draws=draws+n;
    if ~ok
        if strcmp(sample_reason,'sample_in_obstacle'), continue; end
        reason=sample_reason; break;
    end
    [nearest,d]=active_nearest(T,sample); origin=T.X(nearest,:);
    q=sample; if d>o.step, q=origin+o.step*(sample-origin)/d; end
    tr=struct('iteration',iter,'growA',growA,'sample',sample, ...
        'nearest_id',nearest,'steer_origin',origin,'new_point',q, ...
        'accepted',false,'new_id',0,'chosen_parent_id',0, ...
        'near_ids',zeros(1,0),'sorted_parent_ids',zeros(1,0), ...
        'near_radius',0,'optimization_active',result.success || ~o.conditional_activation, ...
        'rewired',zeros(0,3),'opposite_nearest_id',0,'connection_center',[], ...
        'connection_candidates',zeros(1,0),'connection_costs',zeros(1,0), ...
        'chosen_connection_id',0,'contracted',zeros(0,3),'prunedA',zeros(1,0),'prunedB',zeros(1,0));
    if o.record_trace && isempty(trace), trace=tr([]); end
    optimize=tr.optimization_active;
    [near,radius]=near_vertices(T,q,nearest,optimize,o,D);
    costs=T.c(near)+vecnorm(T.X(near,:)-q,2,2);
    [~,order]=sortrows([costs(:) near(:)],[1 2]); ordered=near(order);
    tr.near_ids=near(:)'; tr.near_radius=radius; tr.sorted_parent_ids=ordered(:)';
    parent=0;
    if d>0
        for v=ordered(:)'
            prospective=T.c(v)+norm(T.X(v,:)-q);
            if prospective+norm(q-target(:)')+1e-12>=best, continue; end
            if ~rr.collision(T.X(v,:),q,s), parent=v; break; end
        end
    end
    if parent~=0
        [T,id]=rr.insert(T,q,parent); T.active(id,1)=true;
        accepted=accepted+1; tr.accepted=true; tr.new_id=id; tr.chosen_parent_id=parent;
        if optimize
            for v=near(:)'
                prospective=T.c(id)+norm(q-T.X(v,:));
                if v~=1 && prospective+1e-12<T.c(v) && ...
                        prospective+norm(T.X(v,:)-target(:)')+1e-12<best && ...
                        ~rr.collision(q,T.X(v,:),s)
                    old=T.p(v); T=rr.reparent(T,v,id); rewirings=rewirings+1;
                    tr.rewired(end+1,:)=[v old id];
                end
            end
        end
        if growA, A=T; else, B=T; end
        [result,best,best_link,first_path,first_iteration]= ...
            observe(A,B,links,result,best,first_path,first_iteration,t0,iter);
        [opp,du]=active_nearest(U,q); center=q;
        if du>o.step, center=U.X(opp,:)+o.step*(q-U.X(opp,:))/du; end
        [others,~]=near_vertices(U,center,opp,optimize,o,D);
        fullcost=U.c(others)+vecnorm(U.X(others,:)-q,2,2)+T.c(id);
        [~,order]=sortrows([fullcost(:) others(:)],[1 2]);
        others=others(order); fullcost=fullcost(order);
        tr.opposite_nearest_id=opp; tr.connection_center=center;
        tr.connection_candidates=others(:)'; tr.connection_costs=fullcost(:)';
        for k=1:numel(others)
            if fullcost(k)+1e-12>=best, break; end
            v=others(k);
            if ~rr.collision(q,U.X(v,:),s)
                if growA, pair=[id v]; else, pair=[v id]; end
                links(end+1,:)=pair; connections=connections+1;
                tr.chosen_connection_id=v; break;
            end
        end
        [result,best,best_link,first_path,first_iteration]= ...
            observe(A,B,links,result,best,first_path,first_iteration,t0,iter);
    end
    if result.success && strcmp(o.stop_mode,'first_solution')
        reason='first_solution';
        if o.record_trace, trace(end+1)=tr; end %#ok<AGROW>
        break;
    end
    if result.success && rand<o.contraction_probability
        [A,B,links,changed,event]=contract(A,B,links,best_link,best,s);
        contractions=contractions+changed; tr.contracted=event;
        [result,best,best_link,first_path,first_iteration]= ...
            observe(A,B,links,result,best,first_path,first_iteration,t0,iter);
    end
    if result.success && o.branch_and_bound
        [A,deadA]=prune(A,s.goal,best); [B,deadB]=prune(B,s.start,best);
        pruned=pruned+numel(deadA)+numel(deadB);
        tr.prunedA=deadA(:)'; tr.prunedB=deadB(:)';
    end
    if o.record_trace, trace(end+1)=tr; end %#ok<AGROW>
end
if result.success
    for k=2:size(result.path,1)
        if rr.collision(result.path(k-1,:),result.path(k,:),s)
            error('jordan:Path','A final path edge collides.');
        end
    end
    assert(isequal(result.path(1,:),s.start(:)') && isequal(result.path(end,:),s.goal(:)'));
end
result.treeA=A; result.treeB=B; result.elapsed=toc(t0);
result.details=struct('source_doi','', ...
    'source_id','MIT-CSAIL-TR-2013-021', ...
    'source_url','https://dspace.mit.edu/entities/publication/8184023e-038f-4cda-94da-32a8c7f2812b', ...
    'options',o,'stop_reason',reason,'trace',trace,'first_path',first_path, ...
    'first_iteration',first_iteration,'final_cost',best,'sampling_draws',draws, ...
    'accepted_samples',accepted,'rewirings',rewirings,'connections',connections, ...
    'contractions',contractions,'pruned_vertices',pruned,'cross_edges',links, ...
    'best_cross_edge',best_link,'active_nodes',[sum(A.active) sum(B.active)], ...
    'native_termination','N iterations; finite time_limit is an explicit runner safeguard.', ...
    'variant','Jordan-Perez sorted connection + conditional activation + bounds + contraction + pruning', ...
    'pseudocode_repairs','See README.md; not a literal author-code replication.');
end

function o=options(s,overrides)
D=numel(s.start); volume=prod(s.space(:,2)-s.space(:,1));
unitball=pi^(D/2)/gamma(1+D/2);
o=struct('max_iter',1000,'time_limit',Inf,'step',Inf,'near_radius',Inf, ...
    'gamma',2*(1+1/D)^(1/D)*(volume/unitball)^(1/D), ...
    'goal_radius',1,'goal_bias',0,'record_trace',false,'stop_mode','native', ...
    'conditional_activation',true,'contraction_probability',.1,'branch_and_bound',true, ...
    'sample_sequence',zeros(0,D),'max_sample_attempts',100000);
fields=fieldnames(overrides);
for k=1:numel(fields)
    if ~isfield(o,fields{k}), error('jordan:Option','Unknown option %s.',fields{k}); end
    o.(fields{k})=overrides.(fields{k});
end
validateattributes(o.max_iter,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(o.time_limit,{'numeric'},{'scalar','positive'});
validateattributes(o.step,{'numeric'},{'scalar','positive'});
validateattributes(o.near_radius,{'numeric'},{'scalar','positive'});
validateattributes(o.gamma,{'numeric'},{'scalar','positive','finite'});
validateattributes(o.goal_bias,{'numeric'},{'scalar','>=',0,'<=',1,'finite'});
validateattributes(o.contraction_probability,{'numeric'},{'scalar','>=',0,'<=',1,'finite'});
validateattributes(o.max_sample_attempts,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(o.conditional_activation,{'numeric','logical'},{'scalar'});
validateattributes(o.branch_and_bound,{'numeric','logical'},{'scalar'});
validateattributes(o.record_trace,{'numeric','logical'},{'scalar'});
if ~ismember(string(o.stop_mode),["native","budget","first_solution"])
    error('jordan:Option','Unknown stop_mode.');
end
if ~isempty(o.sample_sequence) && (size(o.sample_sequence,2)~=D || ...
        any(~isfinite(o.sample_sequence),'all'))
    error('jordan:Option','sample_sequence must contain finite D-dimensional rows.');
end
end

function [q,ok,reason,draws]=draw_sample(s,o,iter,target,t0)
ok=false; reason=''; draws=0; q=zeros(1,numel(s.start));
if ~isempty(o.sample_sequence)
    if iter>size(o.sample_sequence,1), reason='sample_sequence_exhausted'; return; end
    q=o.sample_sequence(iter,:); draws=1;
    if rr.collision(q,q,s), reason='sample_in_obstacle'; else, ok=true; end
    return;
end
for draws=1:o.max_sample_attempts
    if toc(t0)>=o.time_limit, reason='time_limit'; return; end
    if rand<o.goal_bias, q=target(:)';
    else, q=s.space(:,1)'+rand(1,numel(s.start)).*(s.space(:,2)-s.space(:,1))'; end
    if ~rr.collision(q,q,s), ok=true; return; end
end
reason='sample_attempt_limit';
end

function [id,d]=active_nearest(T,q)
ids=find(T.active); [d,k]=min(vecnorm(T.X(ids,:)-q,2,2)); id=ids(k);
end

function [ids,radius]=near_vertices(T,q,nearest,optimize,o,D)
if ~optimize, ids=nearest; radius=0; return; end
n=max(2,sum(T.active)); radius=min(o.near_radius,o.gamma*(log(n)/n)^(1/D));
ids=find(T.active & vecnorm(T.X-q,2,2)<=radius);
end

function [r,best,best_link,first_path,first_iter]=observe(A,B,links,r,best,first_path,first_iter,t0,iter)
best_link=0;
if isempty(links), return; end
valid=A.active(links(:,1)) & B.active(links(:,2)); indices=find(valid);
if isempty(indices), return; end
pairs=links(indices,:);
costs=A.c(pairs(:,1))+B.c(pairs(:,2))+vecnorm(A.X(pairs(:,1),:)-B.X(pairs(:,2),:),2,2);
[cost,k]=min(costs); best_link=indices(k);
if cost+1e-12<best
    pair=links(best_link,:); P=rr.path(A,pair(1)); Q=flipud(rr.path(B,pair(2)));
    if isequal(P(end,:),Q(1,:)), Q(1,:)=[]; end
    P=[P;Q]; now=toc(t0);
    if ~r.success, r.first_time=now; first_path=P; first_iter=iter; end
    r.success=true; r.path=P; r.raw_path=P; best=cost;
    r.history(end+1,:)=[iter now cost];
end
end

function ids=chain(T,id)
ids=zeros(1,0);
while id~=0, ids(end+1)=id; id=T.p(id); end %#ok<AGROW>
ids=fliplr(ids);
end

function [A,B,links,changed,event]=contract(A,B,links,k,best,s)
changed=0; event=zeros(0,3); if k==0, return; end
ai=chain(A,links(k,1)); bi=fliplr(chain(B,links(k,2)));
tags=[ones(numel(ai),1) ai(:); 2*ones(numel(bi),1) bi(:)];
points=[A.X(ai,:);B.X(bi,:)];
if size(points,1)<3, return; end
p=sort(randperm(size(points,1),2)); i=p(1); j=p(2);
if j==i+1 || rr.collision(points(i,:),points(j,:),s), return; end
if tags(i,1)==tags(j,1)
    if tags(i,1)==1
        parent=tags(i,2); child=tags(j,2);
        if A.c(parent)+norm(points(i,:)-points(j,:))+1e-12>=A.c(child), return; end
        A=rr.reparent(A,child,parent); event=[1 child parent];
    else
        parent=tags(j,2); child=tags(i,2);
        if B.c(parent)+norm(points(i,:)-points(j,:))+1e-12>=B.c(child), return; end
        B=rr.reparent(B,child,parent); event=[2 child parent];
    end
else
    a=tags(i,2); b=tags(j,2);
    cost=A.c(a)+B.c(b)+norm(points(i,:)-points(j,:));
    if cost+1e-12>=best, return; end
    links(end+1,:)=[a b]; event=[3 a b];
end
changed=1;
end

function [T,dead]=prune(T,target,best)
bad=T.active & T.c+vecnorm(T.X-target(:)',2,2)>best+1e-10;
bad(1)=false;
while true
    next=bad;
    for id=2:numel(T.p), if bad(T.p(id)), next(id)=true; end; end
    if isequal(next,bad), break; end
    bad=next;
end
dead=find(T.active & bad); T.active(dead)=false;
end
