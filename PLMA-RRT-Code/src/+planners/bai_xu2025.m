function out = bai_xu2025(s,options)
if nargin<2, options=struct; end
D=numel(s.start); s.start=s.start(:)'; s.goal=s.goal(:)';
if ~ismember(D,[2 3]) || numel(s.goal)~=D
    error('bai:Dimension','BAI supports matching 2D or 3D endpoints.');
end
opt=configuration(s,options); timer=tic;
if rr.collision(s.start,s.start,s) || rr.collision(s.goal,s.goal,s)
    error('bai:Endpoint','Start and goal must be collision free.');
end
A=rr.new_tree(s.start); B=rr.new_tree(s.goal);
best=zeros(0,D); cbest=Inf; first=NaN; history=zeros(0,3); bridges=zeros(0,2);
trace=struct('iteration',{},'phase',{},'region',{},'theta',{},'d1',{}, ...
    'd2',{},'sample',{},'candidates',{},'accepted',{},'rewired',{},'best_cost',{});
stats=struct('uniform_samples',0,'rectangle_samples',0,'ellipsoid_samples',0, ...
    'goal_samples',0,'proposals',0,'accepted',0,'angle_rejections',0, ...
    'collision_rejections',0,'rewires',0,'connection_tests',0);
it=0; reason='max_iter';
if norm(s.start-s.goal)<=opt.numeric_tolerance
    best=s.start; cbest=0; first=toc(timer); history=[0 first 0]; reason='coincident_endpoints';
else
    for iteration=1:opt.max_iter
        if toc(timer)>=opt.time_limit, reason='time_limit'; break; end
        it=iteration;
        optimizing=~isempty(best);
        [q,region,theta,d1,d2]=sample(s,best,opt,mod(it,2)==0);
        stats.([region '_samples'])=stats.([region '_samples'])+1;
        ia=rr.nearest(A,q); ib=rr.nearest(B,q);
        qa=A.X(ia,:); qb=B.X(ib,:);
        proposals=[steer(qa,q,opt.step);steer(qb,q,opt.step); ...
            steer(qa,qb,opt.step);steer(qb,qa,opt.step)];
        idsA=[]; idsB=[]; accepted=0; rewired=0;
        order=[1 2 3 4]; if mod(it,2)==0, order=[2 1 4 3]; end
        for k=order
            if toc(timer)>=opt.time_limit, break; end
            stats.proposals=stats.proposals+1;
            if mod(k,2)==1
                [A,id,status,nr]=extend(A,proposals(k,:),ia,s,opt,optimizing,opt.start_heading);
                if id>0, idsA(end+1)=id; end %#ok<AGROW>
            else
                [B,id,status,nr]=extend(B,proposals(k,:),ib,s,opt,optimizing,opt.goal_heading);
                if id>0, idsB(end+1)=id; end %#ok<AGROW>
            end
            accepted=accepted+(id>0); rewired=rewired+nr;
            if strcmp(status,'angle'), stats.angle_rejections=stats.angle_rejections+1; end
            if strcmp(status,'collision'), stats.collision_rejections=stats.collision_rejections+1; end
        end
        stats.accepted=stats.accepted+accepted; stats.rewires=stats.rewires+rewired;
        ca=unique(idsA); cb=unique(idsB);
        if any(sqrt(sum((B.X(idsB,:)-qa).^2,2))<=opt.numeric_tolerance), ca=unique([ca ia]); end
        if any(sqrt(sum((A.X(idsA,:)-qb).^2,2))<=opt.numeric_tolerance), cb=unique([cb ib]); end
        for a=ca(:)'
            for b=cb(:)'
                if ~ismember(a,idsA) && ~ismember(b,idsB), continue; end
                if norm(A.X(a,:)-B.X(b,:))>opt.goal_radius, continue; end
                stats.connection_tests=stats.connection_tests+1;
                if can_join(A,a,B,b,s,opt)
                    bridges(end+1,:)=[a b]; %#ok<AGROW>
                end
            end
        end
        if ~isempty(bridges)
            bridges=unique(bridges,'rows','stable');
            edge=sqrt(sum((A.X(bridges(:,1),:)-B.X(bridges(:,2),:)).^2,2));
            costs=A.c(bridges(:,1))+B.c(bridges(:,2))+edge;
            candidates=find(costs<cbest-opt.numeric_tolerance);
            [~,rank]=sort(costs(candidates));
            for j=candidates(rank)'
                a=bridges(j,1); b=bridges(j,2);
                if ~can_join(A,a,B,b,s,opt), continue; end
                P=deduplicate([rr.path(A,a);flipud(rr.path(B,b))],opt.numeric_tolerance);
                if ~valid_path(P,s,opt), continue; end
                best=P; cbest=path_length(P);
                now=toc(timer); if isnan(first), first=now; end
                history(end+1,:)=[it now cbest]; %#ok<AGROW>
                break;
            end
        end
        if opt.record_trace
            if optimizing, phase='refinement'; else, phase='initial'; end
            trace(end+1)=struct('iteration',it,'phase',phase,'region',region, ...
                'theta',theta,'d1',d1,'d2',d2,'sample',q,'candidates',proposals, ...
                'accepted',accepted,'rewired',rewired,'best_cost',cbest); %#ok<AGROW>
        end
        if ~isempty(best) && strcmp(opt.stop_mode,'first_solution')
            reason='first_solution'; break;
        end
    end
end
search_elapsed=toc(timer); postprocess_elapsed=0;
processed=best; post=struct('status','no_path','sign_correction',true);
if ~isempty(best)
    if opt.postprocess
        postopt=opt; postopt.time_limit=opt.postprocess_time_limit;
        posttimer=tic;
        [processed,post]=postprocess(best,s,postopt,posttimer);
        postprocess_elapsed=toc(posttimer);
    else
        post=struct('status','disabled','sign_correction',true);
    end
end
out=struct('path',processed,'raw_path',best,'processed_path',processed, ...
    'treeA',A,'treeB',B,'iterations',it,'success',~isempty(best), ...
    'first_time',first,'elapsed',toc(timer),'history',history);
out.details=struct('source_doi','10.1016/j.oceaneng.2025.121511', ...
    'options',opt,'trace',trace,'stop_reason',reason,'statistics',stats, ...
    'search_elapsed',search_elapsed,'postprocess_elapsed',postprocess_elapsed, ...
    'time_limit_scope','Search only; native postprocessing uses an independent budget; elapsed includes both', ...
    'postprocess',post, ...
    'raw_length',path_length_or_inf(best),'processed_length',path_length_or_inf(processed), ...
    'dimension',D,'paper_dimension',2,'is_3d_generalization',D==3, ...
    'source_status','Full PDF methods checked; disclosed implementation choices; not original ENC experiment', ...
    'gradient_update','X_next = X - gamma * gradient(L); corrects printed Eq.10 sign', ...
    'refinement_choice','Continue both trees with minimum-cost parent selection and angle-preserving RRT* rewiring', ...
    'sampling_threshold_origin','0.309495 from 2D Eq.7; retained as heuristic in 3D', ...
    'parameter_provenance',parameter_provenance());
end

function o=configuration(s,input)
scale=norm(s.space(:,2)-s.space(:,1));
o=struct('max_iter',2500,'time_limit',Inf,'step',scale/25, ...
    'near_radius',[],'goal_radius',[],'goal_bias',0,'record_trace',false, ...
    'stop_mode','native','theta_max',35*pi/180,'adaptive_threshold',0.309495, ...
    'start_heading',[],'goal_heading',[],'postprocess',true,'smooth',true, ...
    'postprocess_time_limit',Inf, ...
    'safety_distance',[],'learning_rate',[],'safety_iterations',100, ...
    'safety_tolerance',[],'interpolation_spacing',[],'boundary_spacing',[], ...
    'boundary_point_cap',4096,'sampling_attempts',100,'backtrack_steps',18, ...
    'numeric_tolerance',max(1,scale)*1e-10);
names=fieldnames(input);
for k=1:numel(names)
    if ~isfield(o,names{k}), error('bai:Option','Unknown BAI option: %s',names{k}); end
    o.(names{k})=input.(names{k});
end
if isempty(o.near_radius), o.near_radius=3*o.step; end
if isempty(o.goal_radius), o.goal_radius=o.step; end
if isempty(o.safety_distance), o.safety_distance=.24*o.step; end
if isempty(o.learning_rate), o.learning_rate=.05*o.safety_distance; end
if isempty(o.safety_tolerance), o.safety_tolerance=1e-6*max(1,scale); end
if isempty(o.interpolation_spacing), o.interpolation_spacing=max(o.safety_distance/2,o.step/20); end
if isempty(o.boundary_spacing), o.boundary_spacing=max(o.safety_distance/2,o.step/20); end
positive={'step','near_radius','goal_radius','interpolation_spacing','boundary_spacing','numeric_tolerance'};
for k=1:numel(positive)
    v=o.(positive{k}); if ~isscalar(v) || ~isfinite(v) || v<=0, error('bai:Option','%s must be positive.',positive{k}); end
end
integers={'max_iter','safety_iterations','boundary_point_cap','sampling_attempts','backtrack_steps'};
for k=1:numel(integers)
    v=o.(integers{k}); if ~isscalar(v) || ~isfinite(v) || v<0 || v~=fix(v), error('bai:Option','%s must be a nonnegative integer.',integers{k}); end
end
if o.boundary_point_cap<16 || o.sampling_attempts<1 || o.backtrack_steps<1 || ...
        ~isscalar(o.time_limit) || isnan(o.time_limit) || o.time_limit<=0 || ...
        ~isscalar(o.postprocess_time_limit) || isnan(o.postprocess_time_limit) || o.postprocess_time_limit<=0 || ...
        ~isscalar(o.theta_max) || ~isfinite(o.theta_max) || o.theta_max<=0 || o.theta_max>pi || ...
        ~isscalar(o.goal_bias) || ~isfinite(o.goal_bias) || o.goal_bias<0 || o.goal_bias>1 || ...
        ~isscalar(o.adaptive_threshold) || ~isfinite(o.adaptive_threshold) || o.adaptive_threshold<=0 || o.adaptive_threshold>=.5 || ...
        ~isscalar(o.safety_distance) || ~isfinite(o.safety_distance) || o.safety_distance<0 || ...
        ~isscalar(o.learning_rate) || ~isfinite(o.learning_rate) || o.learning_rate<0 || ...
        ~isscalar(o.safety_tolerance) || ~isfinite(o.safety_tolerance) || o.safety_tolerance<0 || ...
        ~ismember(o.stop_mode,{'native','first_solution'})
    error('bai:Option','Invalid BAI bounds, counts, stopping mode, angle or safety option.');
end
for name={'start_heading','goal_heading'}
    v=o.(name{1});
    if ~isempty(v)
        if numel(v)~=numel(s.start) || any(~isfinite(v)) || norm(v)==0
            error('bai:Option','Headings must be nonzero D-dimensional outward tree direction vectors.');
        end
        o.(name{1})=v(:)'/norm(v);
    end
end
end

function p=parameter_provenance()
p=struct('theta_max','35 degrees, Table 2 p10 and Table 4 p12', ...
    'max_iter','2500, Table 4 p12; unified override allowed', ...
    'adaptive_threshold','0.309495, Eq.7 p6', ...
    'step','Default workspace diagonal/25 is geometric adaptation; paper uses passage width 250/270/450 m', ...
    'safety_distance','Default .24*step uses Table4 archipelago ratio 60/250; not a universal vessel constant', ...
    'goal_bias','0 by default (uniform pre-solution sampling); nonzero is explicit user override', ...
    'unspecified','Near/connection radii, gamma, GD stop, boundary resolution, spline knots and generic root headings are not fully specified by paper');
end

function [q,kind,theta,d1,d2]=sample(s,P,o,reverse)
D=numel(s.start); theta=NaN; d1=NaN; d2=NaN;
if o.goal_bias>0 && rand<o.goal_bias
    kind='goal'; if reverse, q=s.start; else, q=s.goal; end; return;
end
if isempty(P)
    kind='uniform'; q=s.space(:,1)'+rand(1,D).*(s.space(:,2)-s.space(:,1))'; return;
end
axis=(s.goal-s.start)'; cmin=norm(axis); axis=axis/cmin;
R=[axis null(axis')]; center=(s.start+s.goal)/2;
coordinates=(P-center)*R; d2=max(abs(coordinates(:,1)));
d1=max(sqrt(sum(coordinates(:,2:end).^2,2))); cost=path_length(P); theta=d1/max(cost,eps);
if theta<=o.adaptive_threshold, kind='rectangle'; else, kind='ellipsoid'; end
for attempt=1:o.sampling_attempts
    if strcmp(kind,'rectangle')
        z=(2*rand(1,D)-1).*[d2 repmat(d1,1,D-1)];
    else
        direction=randn(1,D); direction=direction/max(norm(direction),eps);
        ball=direction*rand^(1/D);
        axes=[cost/2 repmat(sqrt(max(0,cost^2-cmin^2))/2,1,D-1)];
        z=ball.*axes;
    end
    q=center+z*R';
    if all(q>=s.space(:,1)') && all(q<=s.space(:,2)'), return; end
end
q=P(randi(size(P,1)),:);
end

function q=steer(a,b,step)
v=b-a; d=norm(v); if d==0, q=a; else, q=a+min(step,d)*v/d; end
end

function [T,id,status,nrewire]=extend(T,q,near,s,o,refine,heading)
id=0; status='duplicate'; nrewire=0;
if norm(q-T.X(near,:))<=o.numeric_tolerance, return; end
if min(sqrt(sum((T.X-q).^2,2)))<=o.numeric_tolerance, return; end
if ~incoming_ok(T,near,q,o,heading), status='angle'; return; end
if rr.collision(T.X(near,:),q,s), status='collision'; return; end
parent=near;
if refine
    distances=sqrt(sum((T.X-q).^2,2)); candidates=find(distances<=o.near_radius);
    [~,idx]=sort(T.c(candidates)+distances(candidates));
    for k=candidates(idx)'
        if T.c(k)+distances(k)>=T.c(parent)+distances(parent), break; end
        if incoming_ok(T,k,q,o,heading) && ~rr.collision(T.X(k,:),q,s)
            parent=k; break;
        end
    end
end
[T,id]=rr.insert(T,q,parent); status='accepted';
if ~refine, return; end
for k=candidates'
    if k==1 || k==parent || T.c(id)+distances(k)>=T.c(k)-o.numeric_tolerance, continue; end
    if ~incoming_ok(T,id,T.X(k,:),o,heading) || rr.collision(q,T.X(k,:),s), continue; end
    kids=find(T.p==k); okay=true;
    for child=kids'
        if ~angle_ok(T.X(k,:)-q,T.X(child,:)-T.X(k,:),o.theta_max), okay=false; break; end
    end
    if okay, T=rr.reparent(T,k,id); nrewire=nrewire+1; end
end
end

function yes=incoming_ok(T,parent,q,o,heading)
if T.p(parent)==0
    yes=isempty(heading) || angle_ok(heading,q-T.X(parent,:),o.theta_max);
else
    yes=angle_ok(T.X(parent,:)-T.X(T.p(parent),:),q-T.X(parent,:),o.theta_max);
end
end

function yes=angle_ok(a,b,limit)
na=norm(a); nb=norm(b);
if na<eps || nb<eps, yes=true; else, yes=dot(a,b)/(na*nb)>=cos(limit)-1e-12; end
end

function yes=can_join(A,a,B,b,s,o)
yes=false; v=B.X(b,:)-A.X(a,:);
if rr.collision(A.X(a,:),B.X(b,:),s), return; end
if norm(v)<=o.numeric_tolerance
    if A.p(a)>0 && B.p(b)>0
        yes=angle_ok(A.X(a,:)-A.X(A.p(a),:),B.X(B.p(b),:)-B.X(b,:),o.theta_max);
    else
        yes=true;
    end
else
    yes=incoming_ok(A,a,B.X(b,:),o,o.start_heading) && incoming_ok(B,b,A.X(a,:),o,o.goal_heading);
end
end

function yes=valid_path(P,s,o)
yes=true;
for k=1:size(P,1)-1
    if rr.collision(P(k,:),P(k+1,:),s), yes=false; return; end
end
for k=2:size(P,1)-1
    if ~angle_ok(P(k,:)-P(k-1,:),P(k+1,:)-P(k,:),o.theta_max), yes=false; return; end
end
if size(P,1)>1
    if ~isempty(o.start_heading) && ~angle_ok(o.start_heading,P(2,:)-P(1,:),o.theta_max), yes=false; end
    if ~isempty(o.goal_heading) && ~angle_ok(o.goal_heading,P(end-1,:)-P(end,:),o.theta_max), yes=false; end
end
end

function [P,report]=postprocess(raw,s,o,timer)
P=raw;
report=struct('status','complete','sign_correction',true,'boundary_points',0, ...
    'loss_before',0,'loss_after',0,'iterations',0,'loss_history',zeros(0,1), ...
    'adjusted_node_count',0,'inserted_node_count',0,'smoothing_attempted',false, ...
    'smoothing_accepted',false,'smoothed_curve_representation','Validated sampled polyline', ...
    'safety_guarantee','Collision-free segments; boundary-point loss reduced; d_ideal is not certified everywhere');
if size(raw,1)<2, report.status='trivial_path'; return; end
if toc(timer)>=o.time_limit, report.status='time_limit_before_postprocess'; return; end
boundary=boundary_points(s,o); report.boundary_points=size(boundary,1);
[P,original]=densify(raw,o.interpolation_spacing); base=P;
report.inserted_node_count=sum(~original);
[loss,grad]=safety_loss(P,boundary,o.safety_distance);
report.loss_before=loss; report.loss_history=loss;
for k=1:o.safety_iterations
    if toc(timer)>=o.time_limit, report.status='time_limit'; break; end
    if loss<=o.safety_tolerance || max(vecnorm(grad,2,2))<=o.safety_tolerance, break; end
    step=o.learning_rate; accepted=false;
    for bt=1:o.backtrack_steps
        proposal=P-step*grad; proposal([1 end],:)=P([1 end],:);
        if valid_path(proposal,s,o)
            [newloss,newgrad]=safety_loss(proposal,boundary,o.safety_distance);
            if newloss<loss-max(1,loss)*1e-12, accepted=true; break; end
        end
        step=step/2;
    end
    if ~accepted, report.status='line_search_stalled'; break; end
    improvement=loss-newloss; P=proposal; loss=newloss; grad=newgrad;
    report.iterations=k; report.loss_history(end+1,1)=loss;
    if improvement<=o.safety_tolerance, break; end
end
report.loss_after=loss;
adjusted=sqrt(sum((P-base).^2,2))>o.numeric_tolerance;
report.adjusted_node_count=sum(adjusted);
keep=original | adjusted; keep([1 end])=true;
for k=find(~keep)'
    before=find(keep(1:k-1),1,'last'); after=k+find(keep(k+1:end),1,'first');
    if isempty(before) || isempty(after) || rr.collision(P(before,:),P(after,:),s), keep(k)=true; end
end
reduced=P(keep,:); moved=adjusted(keep);
if valid_path(reduced,s,o), P=reduced; else, moved=adjusted; end
if o.smooth && size(P,1)>=3 && toc(timer)<o.time_limit
    report.smoothing_attempted=true;
    u=linspace(0,1,size(P,1)); guide=bspline_guide(P,u);
    anchors=moved; anchors([1 end])=true; guide(anchors,:)=P(anchors,:);
    count=max(size(P,1),ceil(path_length(P)/o.interpolation_spacing)*3+1);
    candidate=spline(u,guide',unique([u linspace(0,1,count)]))';
    candidate([1 end],:)=P([1 end],:);
    checkbefore=arc_resample(P,max(40,count)); checkafter=arc_resample(candidate,max(40,count));
    lb=safety_loss(checkbefore,boundary,o.safety_distance);
    la=safety_loss(checkafter,boundary,o.safety_distance);
    if valid_path(candidate,s,o) && la<=lb+o.safety_tolerance
        P=candidate; report.smoothing_accepted=true;
    end
end
if ~valid_path(P,s,o), error('bai:Postprocess','Postprocessing lost a path invariant.'); end
end

function [P,original]=densify(raw,spacing)
P=raw(1,:); original=true;
for k=1:size(raw,1)-1
    count=max(1,ceil(norm(raw(k+1,:)-raw(k,:))/spacing));
    P=[P;raw(k,:)+(1:count)'/count.*(raw(k+1,:)-raw(k,:))]; %#ok<AGROW>
    original=[original;false(count-1,1);true]; %#ok<AGROW>
end
end

function [loss,grad]=safety_loss(P,boundary,safe)
grad=zeros(size(P)); loss=0;
if isempty(boundary) || safe==0, return; end
for i=1:size(P,1)
    differences=P(i,:)-boundary; distances=sqrt(sum(differences.^2,2));
    active=distances<safe; loss=loss+sum(safe-distances(active));
    if i>1 && i<size(P,1) && any(active)
        grad(i,:)=-sum(differences(active,:)./max(distances(active),eps),1);
    end
end
end

function B=boundary_points(s,o)
D=size(s.space,1); B=zeros(0,D); h=o.boundary_spacing; cap=o.boundary_point_cap;
for k=1:numel(s.obstacles)
    ob=s.obstacles{k}; c=ob(1:D); points=zeros(0,D);
    if numel(ob)==D+1
        r=ob(end);
        if D==2
            n=min(cap,max(24,ceil(2*pi*r/h))); t=(0:n-1)'*(2*pi/n);
            points=c+r*[cos(t) sin(t)];
        else
            n=min(cap,max(64,ceil(4*pi*r^2/h^2))); j=(0:n-1)';
            z=1-2*(j+.5)/n; t=pi*(3-sqrt(5))*j;
            points=c+r*[sqrt(1-z.^2).*cos(t) sqrt(1-z.^2).*sin(t) z];
        end
    elseif numel(ob)==2*D
        widths=ob(D+1:end); counts=max(2,ceil(widths/h)+1);
        maxSide=max(2,floor((cap/(2*D))^(1/(D-1)))); counts=min(counts,maxSide);
        for axis=1:D
            other=setdiff(1:D,axis); values=cell(1,D-1);
            for j=1:D-1, values{j}=linspace(-widths(other(j))/2,widths(other(j))/2,counts(other(j))); end
            if D==2, coordinates=values{1}(:); else, [g1,g2]=ndgrid(values{:}); coordinates=[g1(:) g2(:)]; end
            for side=[-1 1]
                face=zeros(size(coordinates,1),D); face(:,axis)=side*widths(axis)/2;
                face(:,other)=coordinates; points=[points;face+c]; %#ok<AGROW>
            end
        end
    elseif D==3 && numel(ob)==5
        r=ob(4); height=ob(5); nt=max(24,ceil(2*pi*r/h)); nz=max(2,ceil(height/h)+1);
        nt=min(nt,max(24,floor(sqrt(cap*2)))); nz=min(nz,max(2,floor(cap/(2*nt))));
        t=(0:nt-1)'*2*pi/nt; z=linspace(-height/2,height/2,nz);
        [ti,zi]=ndgrid(t,z); points=c+[r*cos(ti(:)) r*sin(ti(:)) zi(:)];
        nr=min(max(2,ceil(r/h)+1),max(2,floor(cap/(4*nt))));
        [ti,ri]=ndgrid(t,linspace(0,r,nr));
        for side=[-1 1], points=[points;c+[ri(:).*cos(ti(:)) ri(:).*sin(ti(:)) repmat(side*height/2,numel(ti),1)]]; end %#ok<AGROW>
    end
    points=unique(points,'rows','stable');
    if size(points,1)>cap, points=points(unique(round(linspace(1,size(points,1),cap))),:); end
    B=[B;points]; %#ok<AGROW>
end
end

function G=bspline_guide(P,u)
n=size(P,1); degree=min(3,n-1);
knots=[zeros(1,degree+1),(1:n-degree-1)/(n-degree),ones(1,degree+1)];
G=zeros(numel(u),size(P,2));
for j=1:numel(u)
    if u(j)>=1, G(j,:)=P(end,:); continue; end
    basis=double(u(j)>=knots(1:end-1) & u(j)<knots(2:end));
    for p=1:degree
        next=zeros(1,numel(basis)-1);
        for k=1:numel(next)
            den1=knots(k+p)-knots(k); den2=knots(k+p+1)-knots(k+1);
            if den1>0, next(k)=next(k)+(u(j)-knots(k))/den1*basis(k); end
            if den2>0, next(k)=next(k)+(knots(k+p+1)-u(j))/den2*basis(k+1); end
        end
        basis=next;
    end
    G(j,:)=basis(1:n)*P;
end
end

function Q=arc_resample(P,count)
P=deduplicate(P,1e-12); t=[0;cumsum(sqrt(sum(diff(P).^2,2)))];
if t(end)==0, Q=repmat(P(1,:),count,1); else, Q=interp1(t,P,linspace(0,t(end),count)); end
end

function P=deduplicate(P,tol)
if size(P,1)>1, P=P([true;sqrt(sum(diff(P).^2,2))>tol],:); end
end

function l=path_length(P)
l=sum(sqrt(sum(diff(P).^2,2)));
end

function l=path_length_or_inf(P)
if isempty(P), l=Inf; else, l=path_length(P); end
end
