function varargout = fq_support(action,varargin)
switch action
    case 'setup', [varargout{1:nargout}]=setup(varargin{:});
    case 'sample', [varargout{1:nargout}]=sample_free(varargin{:});
    case 'observe', [varargout{1:nargout}]=observe(varargin{:});
    case 'finish', [varargout{1:nargout}]=finish(varargin{:});
    otherwise, error('fq:Action','Unknown F/Q helper action.');
end
end

function [o,T,r,state,t0]=setup(s,overrides,kind)
D=numel(s.start); t0=tic;
o=struct('max_iter',1000,'time_limit',Inf,'step',Inf,'near_radius',10, ...
    'goal_radius',1,'goal_bias',0,'record_trace',false,'stop_mode','native', ...
    'sample_sequence',zeros(0,D),'max_sample_attempts',100000);
if strcmp(kind,'F'), o.dichotomy_distance=2; else, o.depth=2; end
names=fieldnames(overrides);
for k=1:numel(names)
    if ~isfield(o,names{k}), error('fq:Option','Unknown option %s.',names{k}); end
    o.(names{k})=overrides.(names{k});
end
validateattributes(o.max_iter,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(o.time_limit,{'numeric'},{'scalar','positive'});
validateattributes(o.step,{'numeric'},{'scalar','positive'});
validateattributes(o.near_radius,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(o.goal_radius,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(o.goal_bias,{'numeric'},{'scalar','>=',0,'<=',1,'finite'});
validateattributes(o.max_sample_attempts,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(o.record_trace,{'numeric','logical'},{'scalar'});
if ~ismember(string(o.stop_mode),["native","first_solution","budget"])
    error('fq:Option','stop_mode must be native, first_solution, or budget.');
end
if strcmp(kind,'F')
    validateattributes(o.dichotomy_distance,{'numeric'},{'scalar','positive','finite'});
else
    validateattributes(o.depth,{'numeric'},{'scalar','integer','nonnegative','finite'});
end
if ~isempty(o.sample_sequence) && (size(o.sample_sequence,2)~=D || ...
        any(~isfinite(o.sample_sequence),'all'))
    error('fq:Option','sample_sequence must contain finite D-dimensional rows.');
end
if rr.collision(s.start,s.start,s) || rr.collision(s.goal,s.goal,s)
    error('fq:Endpoints','Start and goal must be in free space.');
end
T=rr.new_tree(s.start);
r=struct('path',zeros(0,D),'raw_path',zeros(0,D),'treeA',T,'treeB',[], ...
    'iterations',0,'success',false,'first_time',NaN,'elapsed',0, ...
    'history',zeros(0,3),'details',struct);
state=struct('goal_ids',zeros(0,1),'goal_edges',zeros(0,1),'best',Inf, ...
    'first_path',zeros(0,D),'first_iteration',NaN);
[r,state]=observe(T,1,s,o,r,state,t0,0);
end

function [q,ok,reason,draws]=sample_free(s,o,iter,t0)
draws=0; ok=false; reason=''; q=zeros(1,numel(s.start));
if ~isempty(o.sample_sequence)
    if iter>size(o.sample_sequence,1), reason='sample_sequence_exhausted'; return; end
    q=o.sample_sequence(iter,:); draws=1;
    if rr.collision(q,q,s), reason='sample_in_obstacle'; else, ok=true; end
    return;
end
for draws=1:o.max_sample_attempts
    if toc(t0)>=o.time_limit, reason='time_limit'; return; end
    if rand<o.goal_bias
        q=s.goal(:)';
    else
        q=s.space(:,1)'+rand(1,numel(s.start)).*(s.space(:,2)-s.space(:,1))';
    end
    if ~rr.collision(q,q,s), ok=true; return; end
end
reason='sample_attempt_limit';
end

function [r,state]=observe(T,new_ids,s,o,r,state,t0,iter)
for id=new_ids(:)'
    d=norm(T.X(id,:)-s.goal(:)');
    if d<=o.goal_radius && ~rr.collision(T.X(id,:),s.goal,s)
        state.goal_ids(end+1,1)=id;
        state.goal_edges(end+1,1)=d;
    end
end
if isempty(state.goal_ids), return; end
[cost,k]=min(T.c(state.goal_ids)+state.goal_edges);
if cost+1e-12<state.best
    P=rr.path(T,state.goal_ids(k));
    if ~isequal(P(end,:),s.goal(:)'), P(end+1,:)=s.goal(:)'; end
    now=toc(t0);
    if ~r.success
        r.first_time=now; state.first_path=P; state.first_iteration=iter;
    end
    state.best=cost; r.success=true; r.path=P; r.raw_path=P;
    r.history(end+1,:)=[iter now cost];
end
end

function r=finish(T,s,o,r,state,t0,trace,details,reason)
if r.success
    if ~isequal(r.path(1,:),s.start(:)') || ~isequal(r.path(end,:),s.goal(:)')
        error('fq:Path','Solution endpoints are inconsistent.');
    end
    for k=2:size(r.path,1)
        if rr.collision(r.path(k-1,:),r.path(k,:),s)
            error('fq:Path','A solution contains a colliding edge.');
        end
    end
end
r.treeA=T; r.elapsed=toc(t0);
details.options=o; details.stop_reason=reason; details.trace=trace;
details.first_path=state.first_path; details.first_iteration=state.first_iteration;
details.final_cost=state.best; details.completed_time=r.elapsed;
details.goal_connection='Certified exact-goal terminal edge within goal_radius; goal excluded from search tree.';
details.steering_adapter=isfinite(o.step);
r.details=details;
end
