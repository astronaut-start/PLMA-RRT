function r=gb_karaman2011(s,options)
if nargin<2, options=struct; end
D=numel(s.start); span=s.space(:,2)-s.space(:,1);
zeta=pi^(D/2)/gamma(D/2+1);
o=struct('max_iter',2500,'time_limit',10,'step',0.04*max(span), ...
 'gamma',2.2*(1+1/D)^(1/D)*(prod(span)/zeta)^(1/D), ...
 'goal_radius',0.05*max(span),'goal_bias',0.1,'record_trace',false, ...
 'stop_mode','native','sample_sequence',zeros(0,D));
for f=fieldnames(options)'
 if ~isfield(o,f{1}), error('gb:Option','Unknown option %s',f{1}); end
 o.(f{1})=options.(f{1});
end
validateattributes(o.step,{'numeric'},{'positive','scalar','finite'});
validateattributes(o.max_iter,{'numeric'},{'nonnegative','integer','scalar'});
validateattributes(o.time_limit,{'numeric'},{'positive','scalar'});
validateattributes(o.goal_bias,{'numeric'},{'scalar','>=',0,'<=',1});
validateattributes(o.gamma,{'numeric'},{'positive','scalar','finite'});
if ~ismember(string(o.stop_mode),["native","budget","first_solution"]), error('gb:Stop','Invalid stop mode'); end
if rr.collision(s.start,s.start,s)||rr.collision(s.goal,s.goal,s), error('gb:Endpoint','Occupied endpoint'); end
t0=tic; T=rr.new_tree(s.start); goals=[]; goalEdges=[]; best=Inf;
P=zeros(0,D); firstP=P; first=NaN; firstN=NaN; history=zeros(0,3);
trace=struct('iteration',{},'sample',{},'nearest_id',{},'nearest',{},'new_point',{},'selected_parent',{},'radius',{},'rewired',{});
reason='max_iter'; it=0;
for iter=1:o.max_iter
 if toc(t0)>=o.time_limit, reason='time_limit'; break; end
 it=iter;
 if ~isempty(o.sample_sequence)
  if iter>size(o.sample_sequence,1), reason='sample_sequence_exhausted'; break; end
  q=o.sample_sequence(iter,:);
 else
  if rand<o.goal_bias
   q=s.goal;
  else
   while true
    q=s.space(:,1)'+rand(1,D).*span';
    if ~rr.collision(q,q,s), break; end
    if toc(t0)>=o.time_limit, break; end
   end
  end
 end
 if rr.collision(q,q,s), continue; end
 [nn,d]=rr.nearest(T,q); origin=T.X(nn,:);
 if d<1e-12, continue; end
 qn=origin+min(1,o.step/d)*(q-origin);
 if rr.collision(origin,qn,s), continue; end
 n=size(T.X,1); rad=min(o.gamma*(log(n)/n)^(1/D),o.step);
 ds=vecnorm(T.X-qn,2,2); near=find(ds<=rad);
 parent=nn; cost=T.c(nn)+norm(qn-origin);
 for j=near'
  c=T.c(j)+ds(j);
  if c<cost && ~rr.collision(T.X(j,:),qn,s), parent=j; cost=c; end
 end
 [T,id]=rr.insert(T,qn,parent); rewired=[];
 for j=near'
  if T.c(id)+ds(j)<T.c(j)-1e-12 && ~rr.collision(qn,T.X(j,:),s)
   T=rr.reparent(T,j,id); rewired(end+1)=j; %#ok<AGROW>
  end
 end
 dg=norm(qn-s.goal);
 if dg<=o.goal_radius && ~rr.collision(qn,s.goal,s), goals(end+1)=id; goalEdges(end+1)=dg; end %#ok<AGROW>
 if ~isempty(goals)
  [c,k]=min(T.c(goals)+goalEdges(:));
  if c<best-1e-12
   best=c; P=rr.path(T,goals(k));
   if norm(P(end,:)-s.goal)>1e-12, P(end+1,:)=s.goal; end
   if isnan(first), first=toc(t0); firstP=P; firstN=size(T.X,1); end
   history(end+1,:)=[iter toc(t0) c]; %#ok<AGROW>
  end
 end
 if o.record_trace
  trace(end+1)=struct('iteration',iter,'sample',q,'nearest_id',nn,'nearest',origin,'new_point',qn,'selected_parent',parent,'radius',rad,'rewired',rewired); %#ok<AGROW>
 end
 if ~isempty(P)&&strcmp(o.stop_mode,'first_solution'), reason='first_solution'; break; end
end
r=struct('path',P,'raw_path',P,'treeA',T,'treeB',[],'iterations',it, ...
 'success',~isempty(P),'first_time',first,'elapsed',toc(t0),'history',history);
r.details=struct('source_doi','10.1177/0278364911406761','options',o,'stop_reason',reason, ...
 'trace',trace,'first_path',firstP,'first_nodes',firstN,'source_status', ...
 'Algorithm 6 + declared fixed goal bias; GB-RRT* is a baseline label, not a uniquely identified paper.');
end
