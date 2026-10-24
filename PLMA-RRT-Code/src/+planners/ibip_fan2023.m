function result = ibip_fan2023(s, overrides)
if nargin<2, overrides=struct; end
o=options(overrides); D=numel(s.start); t0=tic;
s.start=s.start(:)'; s.goal=s.goal(:)';
if numel(s.goal)~=D || rr.collision(s.start,s.start,s) || rr.collision(s.goal,s.goal,s)
    error('fan:Endpoints','Start and goal must be in free space with equal dimensions.');
end
TA=rr.new_tree(s.start); TB=rr.new_tree(s.goal);
result=struct('path',zeros(0,D),'raw_path',zeros(0,D), ...
    'processed_path',zeros(0,D),'treeA',TA,'treeB',TB,'iterations',0, ...
    'success',false,'first_time',NaN,'elapsed',0,'history',zeros(0,3), ...
    'details',struct);
trace=struct([]); bridges=zeros(0,2); best_cost=Inf;
accepted=0; rewirings=0; proposals=0; collisions=0; draws=0;
stop_reason='iteration_limit';
if isequal(s.start,s.goal)
    result.success=true; result.raw_path=s.start; result.first_time=toc(t0);
    result.history=[0,result.first_time,0]; best_cost=0;
    stop_reason='start_equals_goal';
end
for iter=1:o.max_iter
    if strcmp(stop_reason,'start_equals_goal'), break; end
    if toc(t0)>=o.time_limit, stop_reason='time_limit'; break; end
    result.iterations=iter;
    growA=mod(iter,2)==1;
    if growA, T=TA; U=TB; target=s.goal; else, T=TB; U=TA; target=s.start; end
    biased=rand<=o.goal_bias;
    if biased
        sample=target;
    else
        [sample,ndraw,timed_out]=free_sample(s,o,t0); draws=draws+ndraw;
        if timed_out, stop_reason='time_limit_in_sampling'; break; end
    end
    [nearest,~]=rr.nearest(T,sample); q=T.X(nearest,:);
    [force,fg,fr,frep,obstacle_id,dobs]=modified_force(q,target,sample,s,o);
    proposals=proposals+1; qnew=q; inserted=false; parent=0;
    fmax=max(abs(force));
    if all(isfinite(force)) && fmax>0
        qnew=q+o.step*force/fmax;
        if ~rr.collision(q,qnew,s)
            [~,duplicate_distance]=rr.nearest(T,qnew);
            if duplicate_distance>64*eps(max([1,abs(qnew)]))
                [T,new_id,changed,used_radius]=star_insert(T,qnew,nearest,s,o);
                parent=T.p(new_id); accepted=accepted+1;
                rewirings=rewirings+changed; inserted=true;
            else
                used_radius=NaN;
            end
        else
            collisions=collisions+1; used_radius=NaN;
        end
    else
        used_radius=NaN;
    end
    if o.record_trace
        row=struct('iteration',iter,'active_tree',1+~growA,'target',target, ...
            'sample',sample,'goal_biased',biased,'nearest_id',nearest, ...
            'nearest',q,'new_point',qnew,'force',force,'goal_attraction',fg, ...
            'sample_attraction',fr,'repulsion',frep,'obstacle_id',obstacle_id, ...
            'obstacle_center_distance',dobs,'step',o.step, ...
            'inserted',inserted,'selected_parent',parent,'near_radius_used',used_radius);
        if isempty(trace), trace=row; else, trace(end+1)=row; end %#ok<AGROW>
    end
    if growA, TA=T; else, TB=T; end
    if ~inserted, continue; end
    [other_id,gap]=rr.nearest(U,qnew);
    if gap<=o.goal_radius
        ds=vecnorm(U.X-qnew,2,2);
        candidates=unique([other_id;find(ds<=o.near_radius)]);
        [~,order]=sort(U.c(candidates)+ds(candidates));
        for b=candidates(order)'
            if ~rr.collision(qnew,U.X(b,:),s)
                if growA, bridges(end+1,:)=[new_id,b]; else, bridges(end+1,:)=[b,new_id]; end %#ok<AGROW>
                break;
            end
        end
    end
    if ~isempty(bridges)
        bc=TA.c(bridges(:,1))+TB.c(bridges(:,2))+ ...
            vecnorm(TA.X(bridges(:,1),:)-TB.X(bridges(:,2),:),2,2);
        [candidate_cost,j]=min(bc);
        if candidate_cost<best_cost-1e-11*max(1,best_cost) || ~result.success
            pa=rr.path(TA,bridges(j,1)); pb=rr.path(TB,bridges(j,2));
            candidate=[pa;flipud(pb)];
            candidate=candidate([true;vecnorm(diff(candidate),2,2)>0],:);
            if path_valid(candidate,s) && toc(t0)<=o.time_limit
                result.raw_path=candidate; best_cost=path_length(candidate);
                if ~result.success, result.first_time=toc(t0); end
                result.success=true;
                result.history(end+1,:)=[iter,toc(t0),best_cost]; %#ok<AGROW>
                if strcmp(o.stop_mode,'first_solution')
                    stop_reason='first_solution'; break;
                end
            end
        end
    end
end
planning_elapsed=toc(t0);
result.path=result.raw_path;
smoothing=struct('attempted',false,'accepted',false,'reason','no_solution', ...
    'continuous_curve_valid',false,'sampled_segments_valid',false, ...
    'raw_length',best_cost,'candidate_length',NaN,'parameterization','chord_length', ...
    'boundary_condition','MATLAB spline not-a-knot','candidate_path',zeros(0,D), ...
    'pp',struct,'seconds',0);
if result.success
    smooth_start=tic;
    [smoothed,smoothing]=smooth_path(result.raw_path,s,o);
    smoothing.seconds=toc(smooth_start);
    result.processed_path=smoothed;
    if smoothing.accepted, result.path=smoothed; end
end
result.treeA=TA; result.treeB=TB; result.elapsed=toc(t0);
result.details=struct('source_doi','10.1016/j.eswa.2022.119137', ...
    'source_algorithm','ImprovedBi-P-RRT*','protocol_label','IBI-P-RRT*', ...
    'reproduction_level','source-guided reconstruction; not author code', ...
    'options',o,'trace',trace,'stop_reason',stop_reason,'proposals',proposals, ...
    'accepted_nodes',accepted,'collision_rejections',collisions, ...
    'rewirings',rewirings,'uniform_draws',draws,'bridges',bridges, ...
    'planning_elapsed',planning_elapsed,'smoothing',smoothing, ...
    'raw_path_length',path_length(result.raw_path), ...
    'output_path_length',path_length(result.path), ...
    'history_metric','raw collision-free path length; excludes spline', ...
    'time_limit_scope','search and first-solution validation; postprocessing reported separately', ...
    'is_3d_extension',D==3);
end

function o=options(user)
o=struct('max_iter',2000,'time_limit',60,'step',2,'near_radius',10, ...
    'goal_radius',2,'goal_bias',.5,'record_trace',false,'stop_mode','native', ...
    'attraction_gain',1,'repulsion_gain',1,'influence_radius',10, ...
    'goal_distance_power',2,'gamma',10,'radius_cap',10, ...
    'spline_samples_per_edge',16);
if ~isstruct(user) || ~isscalar(user), error('fan:Options','Options must be one struct.'); end
names=fieldnames(user);
for k=1:numel(names)
    if ~isfield(o,names{k}), error('fan:Options','Unknown option: %s',names{k}); end
    o.(names{k})=user.(names{k});
end
if ~isfield(user,'goal_radius'), o.goal_radius=o.step; end
if ~isfield(user,'gamma'), o.gamma=o.near_radius; end
if ~isfield(user,'radius_cap'), o.radius_cap=o.near_radius; end
positive={'step','near_radius','goal_radius','influence_radius','gamma','radius_cap'};
for k=1:numel(positive)
    v=o.(positive{k});
    if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || v<=0
        error('fan:Options','%s must be finite and positive.',positive{k});
    end
end
for key={'attraction_gain','repulsion_gain'}
    v=o.(key{1});
    if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || v<0
        error('fan:Options','%s must be finite and nonnegative.',key{1});
    end
end
for key={'max_iter','goal_distance_power','spline_samples_per_edge'}
    v=o.(key{1});
    if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || v~=fix(v) || v<1
        error('fan:Options','%s must be a positive integer.',key{1});
    end
end
if ~isscalar(o.time_limit) || isnan(o.time_limit) || o.time_limit<=0 || ...
        ~isscalar(o.goal_bias) || ~isfinite(o.goal_bias) || o.goal_bias<0 || o.goal_bias>1
    error('fan:Options','Invalid time limit or goal bias.');
end
if ~(islogical(o.record_trace) && isscalar(o.record_trace))
    error('fan:Options','record_trace must be a logical scalar.');
end
if ~any(strcmp(o.stop_mode,{'native','first_solution'}))
    error('fan:Options','stop_mode must be native or first_solution.');
end
end

function [x,draws,timed_out]=free_sample(s,o,t0)
draws=0; timed_out=false; x=zeros(1,size(s.space,1));
while true
    if toc(t0)>=o.time_limit, timed_out=true; return; end
    draws=draws+1; x=s.space(:,1)'+rand(1,size(s.space,1)).*diff(s.space,1,2)';
    if ~rr.collision(x,x,s), return; end
end
end

function [F,Fg,Fr,Frep,idx,dobs]=modified_force(q,target,sample,s,o)
Fg=o.attraction_gain*(target-q); Fr=o.attraction_gain*(sample-q);
Frep=zeros(size(q)); idx=0; dobs=Inf;
for k=1:numel(s.obstacles)
    obstacle=s.obstacles{k}; delta=q-obstacle(1:numel(q)); d=norm(delta);
    if d<dobs, dobs=d; idx=k; end
end
if idx>0 && dobs<=o.influence_radius && dobs>0
    obs=s.obstacles{idx}; unitOR=(q-obs(1:numel(q)))/dobs;
    dg=norm(target-q); unitRG=zeros(size(q));
    if dg>0, unitRG=(target-q)/dg; end
    c=1/dobs-1/o.influence_radius; n=o.goal_distance_power;
    F1=o.repulsion_gain*c*dg^n/dobs^2;
    F2=(n/2)*o.repulsion_gain*c^2*dg^(n-1);
    Frep=F1*unitOR+F2*unitRG;
end
F=Fg+Fr+Frep;
end

function [T,new_id,changed,radius]=star_insert(T,q,nearest,s,o)
n=size(T.X,1); D=size(T.X,2);
radius=min(o.gamma*(log(max(n,2))/max(n,2))^(1/D),o.radius_cap);
ds=vecnorm(T.X-q,2,2); near=find(ds<=radius);
if isempty(near), radius=o.near_radius; near=find(ds<=radius); end
parent=nearest; cost=T.c(nearest)+ds(nearest);
for k=near'
    proposed=T.c(k)+ds(k);
    if proposed<cost && ~rr.collision(T.X(k,:),q,s), parent=k; cost=proposed; end
end
[T,new_id]=rr.insert(T,q,parent); changed=0;
for k=near'
    if T.p(k)==0 || k==parent, continue; end
    if T.c(new_id)+ds(k)<T.c(k)-1e-11 && ~rr.collision(q,T.X(k,:),s)
        T=rr.reparent(T,k,new_id); changed=changed+1;
    end
end
end

function good=path_valid(P,s)
good=~isempty(P) && isequal(P(1,:),s.start) && isequal(P(end,:),s.goal);
if ~good, return; end
for k=1:max(1,size(P,1)-1)
    if rr.collision(P(k,:),P(min(k+1,size(P,1)),:),s), good=false; return; end
end
end

function L=path_length(P)
if isempty(P), L=Inf; else, L=sum(vecnorm(diff(P),2,2)); end
end

function [P,report]=smooth_path(raw,s,o)
D=size(raw,2); P=zeros(0,D);
report=struct('attempted',true,'accepted',false,'reason','', ...
    'continuous_curve_valid',false,'sampled_segments_valid',false, ...
    'raw_length',path_length(raw),'candidate_length',NaN, ...
    'parameterization','chord_length','boundary_condition','MATLAB spline not-a-knot', ...
    'candidate_path',zeros(0,D),'pp',struct,'seconds',0);
if size(raw,1)<=2
    P=raw; report.accepted=true; report.reason='linear_or_single_point';
    report.continuous_curve_valid=true; report.sampled_segments_valid=true;
    report.candidate_path=P; report.candidate_length=path_length(P); return;
end
t=[0;cumsum(vecnorm(diff(raw),2,2))];
pp=spline(t,raw'); report.pp=pp;
tt=[];
for k=1:numel(t)-1
    ts=linspace(t(k),t(k+1),o.spline_samples_per_edge+1);
    tt=[tt,ts(1:end-1)]; %#ok<AGROW>
end
tt=[tt,t(end)]; candidate=ppval(pp,tt)';
candidate(1,:)=raw(1,:); candidate(end,:)=raw(end,:);
report.candidate_path=candidate; report.candidate_length=path_length(candidate);
report.continuous_curve_valid=curve_valid(pp,s);
report.sampled_segments_valid=path_valid(candidate,s);
if report.continuous_curve_valid && report.sampled_segments_valid
    P=candidate; report.accepted=true; report.reason='validated_cubic';
else
    report.reason='spline_collision_or_workspace_excursion_raw_retained';
end
end

function valid=curve_valid(pp,s)
valid=true; D=pp.dim;
for k=1:pp.pieces
    h=pp.breaks(k+1)-pp.breaks(k);
    C=pp.coefs((k-1)*D+(1:D),:);
    C=C.*(h.^(size(C,2)-1:-1:0));
    for dim=1:D
        cuts=[0;1;unit_roots(polyder(C(dim,:)))];
        v=polyval(C(dim,:),cuts);
        if min(v)<s.space(dim,1) || max(v)>s.space(dim,2)
            valid=false; return;
        end
    end
    for j=1:numel(s.obstacles)
        obs=s.obstacles{j}; C0=C; C0(:,end)=C0(:,end)-obs(1:D)';
        if numel(obs)==D+1
            radial=squared_norm_polynomial(C0);
            candidates=[0;1;unit_roots(polyder(radial))];
            d2=polyval(radial,candidates); r2=obs(D+1)^2;
            tol=1e-9*max([1,r2,max(abs(radial))]);
            if min(d2)<=r2+tol, valid=false; return; end
        elseif numel(obs)==2*D
            cuts=[0;1]; half=obs(D+1:2*D)/2;
            for dim=1:D
                lo=C0(dim,:); lo(end)=lo(end)+half(dim);
                hi=C0(dim,:); hi(end)=hi(end)-half(dim);
                cuts=[cuts;unit_roots(lo);unit_roots(hi);unit_roots(polyder(C0(dim,:)))]; %#ok<AGROW>
            end
            cuts=partition_points(cuts); points=zeros(numel(cuts),D);
            for dim=1:D, points(:,dim)=polyval(C0(dim,:),cuts); end
            tol=1e-9*max([1,abs(obs)]);
            if any(all(abs(points)<=half+tol,2)), valid=false; return; end
        elseif D==3 && numel(obs)==5
            radial=squared_norm_polynomial(C0(1:2,:)); radial(end)=radial(end)-obs(4)^2;
            lo=C0(3,:); lo(end)=lo(end)+obs(5)/2;
            hi=C0(3,:); hi(end)=hi(end)-obs(5)/2;
            cuts=partition_points([0;1;unit_roots(radial);unit_roots(lo);unit_roots(hi); ...
                unit_roots(polyder(radial));unit_roots(polyder(C0(3,:)))]);
            rv=polyval(radial,cuts); zv=polyval(C0(3,:),cuts);
            tol=1e-9*max([1,abs(obs),abs(radial)]);
            if any(rv<=tol & abs(zv)<=obs(5)/2+tol), valid=false; return; end
        else
            error('fan:Obstacle','Unsupported obstacle representation.');
        end
    end
end
end

function p=squared_norm_polynomial(C)
p=zeros(1,2*size(C,2)-1);
for j=1:size(C,1), p=p+conv(C(j,:),C(j,:)); end
end

function r=unit_roots(p)
first=find(p~=0,1,'first');
if isempty(first), r=zeros(0,1); return; end
p=p(first:end);
if numel(p)<=1, r=zeros(0,1); return; end
r=roots(p); r=real(r(abs(imag(r))<=1e-7 & real(r)>=-1e-10 & real(r)<=1+1e-10));
r=min(1,max(0,r(:)));
end

function p=partition_points(cuts)
cuts=unique(sort(cuts)); p=[cuts;(cuts(1:end-1)+cuts(2:end))/2];
end
