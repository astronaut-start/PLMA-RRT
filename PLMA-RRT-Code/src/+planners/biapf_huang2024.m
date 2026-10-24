function result = biapf_huang2024(s, overrides)
if nargin<2, overrides=struct; end
o=planners.defaults_huang2024(s,overrides);
t0=tic; D=numel(s.start); s.start=s.start(:)'; s.goal=s.goal(:)';
if rr.collision(s.start,s.start,s) || rr.collision(s.goal,s.goal,s)
    error('huang:Endpoints','Start and goal must be in free space.');
end
TA=rr.new_tree(s.start); TB=rr.new_tree(s.goal);
result=struct('path',zeros(0,D),'raw_path',zeros(0,D), ...
    'processed_path',zeros(0,D),'treeA',TA,'treeB',TB,'iterations',0, ...
    'success',false,'first_time',NaN,'elapsed',0,'history',zeros(0,3), ...
    'details',struct);
trace=struct([]); proposals=0; accepted=0; collisions=0; switches=0;
sampling_draws=0; rewirings=0; stop_reason='iteration_limit';
growA=true;
if isequal(s.start,s.goal)
    result.path=s.start; result.raw_path=s.start; result.success=true;
    result.first_time=toc(t0); result.history=[0,result.first_time,0];
    stop_reason='start_equals_goal';
end
for iter=1:o.max_iter
    if result.success, break; end
    if toc(t0)>=o.time_limit, stop_reason='time_limit'; break; end
    result.iterations=iter;
    if growA, T=TA; U=TB; nominal_target=s.goal; else, T=TB; U=TA; nominal_target=s.start; end
    switched=rand<o.target_switch_probability;
    if switched
        target=uniform_sample(s.space); switches=switches+1;
    else
        target=nominal_target;
    end
    current=T.X(end,:);
    pgoal=o.goal_bias_initial*exp(-o.goal_bias_distance_decay*norm(current-target)) ...
        *exp(-o.goal_bias_iteration_decay*iter);
    goal_biased=rand<pgoal;
    if goal_biased
        sample=target;
    else
        [sample,draws,timed_out]=regional_sample(s,o,t0);
        sampling_draws=sampling_draws+draws;
        if timed_out, stop_reason='time_limit_in_sampling'; break; end
    end
    [nearest,~]=rr.nearest(T,sample); q=T.X(nearest,:);
    [F,Nobs,density,Kadj,Fatt,Frep]=improved_force(q,target,s.obstacles,iter,o);
    h=o.step_max*exp(-o.step_density_decay*Nobs);
    proposals=proposals+1; qnew=q; inserted=false; parent=0;
    if all(isfinite(F)) && norm(F)>0 && h>0
        qnew=q+h*F/norm(F);
        if ~rr.collision(q,qnew,s)
            [~,nearest_distance]=rr.nearest(T,qnew);
            if nearest_distance>64*eps(max([1,abs(qnew)]))
                [T,new_id,rewired]=star_insert(T,qnew,nearest,s,o);
                inserted=true; parent=T.p(new_id);
                accepted=accepted+1; rewirings=rewirings+rewired;
            end
        else
            collisions=collisions+1;
        end
    end
    if o.record_trace
        row=struct('iteration',iter,'active_tree',1+~growA, ...
            'target',target,'sample',sample,'nearest_id',nearest,'nearest',q, ...
            'new_point',qnew,'force',F,'attraction',Fatt,'base_repulsion',Frep, ...
            'bias_reference_point',current,'nearby_obstacles',Nobs,'density',density, ...
            'repulsion_coefficient',Kadj,'step',h,'p_goal',pgoal, ...
            'target_switched',switched,'goal_biased',goal_biased, ...
            'inserted',inserted,'selected_parent',parent);
        if isempty(trace), trace=row; else, trace(end+1)=row; end %#ok<AGROW>
    end
    if growA, TA=T; else, TB=T; end
    if ~inserted, continue; end
    [other_id,bridge_length]=rr.nearest(U,qnew);
    if bridge_length<=o.connect_radius && ~rr.collision(qnew,U.X(other_id,:),s)
        if growA
            pa=rr.path(T,new_id); pb=rr.path(U,other_id);
        else
            pa=rr.path(U,other_id); pb=rr.path(T,new_id);
        end
        candidate=[pa;flipud(pb)];
        if size(candidate,1)>1
            candidate=candidate([true;vecnorm(diff(candidate),2,2)>0],:);
        end
        valid=isequal(candidate(1,:),s.start(:)') && isequal(candidate(end,:),s.goal(:)');
        for edge=1:size(candidate,1)-1
            if rr.collision(candidate(edge,:),candidate(edge+1,:),s), valid=false; break; end
        end
        if valid && toc(t0)<=o.time_limit
            result.path=candidate; result.raw_path=candidate;
            result.success=true; result.first_time=toc(t0);
            result.history=[iter,result.first_time,sum(vecnorm(diff(candidate),2,2))];
            stop_reason='first_valid_connection'; break;
        end
    end
    growA=~growA;
end
result.treeA=TA; result.treeB=TB; result.elapsed=toc(t0);
result.details=struct('source_doi','10.3390/drones8120760','source_version','v2, 26 pages', ...
    'status','source-guided reconstruction; unpublished parameters explicitly chosen', ...
    'native_stop','first successful connection (Figure 3)', ...
    'stop_reason',stop_reason,'options',o,'proposals',proposals,'accepted',accepted, ...
    'collision_rejections',collisions,'switches',switches,'sampling_draws',sampling_draws, ...
    'rewirings',rewirings,'trace',trace, ...
    'connection_rule','nearest passive-tree node; full segment checked', ...
    'finite_connection_override',isfinite(o.connect_radius), ...
    'alternation_rule','swap only after an accepted collision-free extension', ...
    'apf_assembly','Figure 3 force expansion direction at nearest node; Algorithm 3 shifted-sample conflict documented', ...
    'postprocess','None specified by Huang; processed_path intentionally empty');
end

function [T,id,rewired] = star_insert(T,q,nearest,s,o)
d=vecnorm(T.X-q,2,2); near=find(d<=o.near_radius);
parent=nearest; best=T.c(nearest)+d(nearest);
for k=near'
    cost=T.c(k)+d(k);
    if cost<best && ~rr.collision(T.X(k,:),q,s), parent=k; best=cost; end
end
[T,id]=rr.insert(T,q,parent); rewired=0;
for k=near'
    if T.p(k)==0 || k==parent, continue; end
    if T.c(id)+d(k)<T.c(k)-1e-12*max(1,T.c(k)) && ~rr.collision(q,T.X(k,:),s)
        T=rr.reparent(T,k,id); rewired=rewired+1;
    end
end
end

function q = uniform_sample(bounds)
q=bounds(:,1)'+rand(1,size(bounds,1)).*(bounds(:,2)-bounds(:,1))';
end

function [sample,draws,timed_out] = regional_sample(s,o,t0)
sample=[]; draws=0; timed_out=false;
a=s.start(:)'; v=s.goal(:)'-a; vv=dot(v,v);
while true
    if toc(t0)>=o.time_limit, timed_out=true; return; end
    q=uniform_sample(s.space); draws=draws+1;
    if vv==0, projected=a; else, projected=a+min(1,max(0,dot(q-a,v)/vv))*v; end
    d=norm(q-projected);
    if rand<exp(-d*d/(2*o.region_sigma^2)), sample=q; return; end
end
end

function [F,Nobs,density,Kadj,Fatt,Frep] = improved_force(q,target,obstacles,iter,o)
Fatt=o.Katt*(target-q); Frep=zeros(size(q)); dgoal=norm(target-q); Nobs=0;
for k=1:numel(obstacles)
    [d,point]=point_to_obstacle(q,obstacles{k});
    if d<o.influence_radius
        Nobs=Nobs+1;
        if d<=0, F=NaN(size(q)); density=NaN; Kadj=NaN; return; end
        factor=1/d-1/o.influence_radius;
        f1=o.Krep*factor*dgoal^o.rep_power/d^2;
        f2=(o.rep_power/2)*o.Krep*factor*dgoal^(o.rep_power-1)/d^2;
        Frep=Frep+(f1+f2)*(q-point)/d;
    end
end
density=Nobs/max(1,numel(obstacles));
if dgoal<o.goal_region_radius
    Kadj=o.Krep*(1-o.density_beta*density)*exp(-0.01*iter);
else
    Kadj=o.Krep*(1+o.density_alpha*density)*exp(-0.01*iter);
end
F=Fatt+Kadj*Frep;
end

function [d,point] = point_to_obstacle(q,o)
D=numel(q); o=o(:)';
if numel(o)==D+1
    v=q-o(1:D); n=norm(v); r=o(end);
    if n<=r, point=q; d=0; else, point=o(1:D)+r*v/n; d=n-r; end
elseif numel(o)==2*D
    point=min(max(q,o(1:D)-o(D+1:end)/2),o(1:D)+o(D+1:end)/2); d=norm(q-point);
elseif D==3 && numel(o)==5
    v=q(1:2)-o(1:2); n=norm(v);
    if n>o(4), xy=o(1:2)+o(4)*v/n; else, xy=q(1:2); end
    point=[xy,min(max(q(3),o(3)-o(5)/2),o(3)+o(5)/2)]; d=norm(q-point);
else
    error('huang:Obstacle','Unsupported obstacle type.');
end
end
