function T = plma_reparent(T, child, parent)
if ~isstruct(T) || isempty(T) || ~isvector(T) || ...
        ~all(isfield(T,{'pos','parent','cost'}))
    error('plma:Tree','T must be a nonempty native PLMA struct vector.');
end
n=numel(T);
if ~valid_index(child,n) || ~valid_index(parent,n)
    error('plma:Parent','Child and parent must index existing nodes.');
end
D=numel(T(1).pos);
if ~ismember(D,[2 3]), error('plma:Position','Positions must be 2D or 3D vectors.'); end
X=zeros(n,D); p=zeros(n,1);
for k=1:n
    q=T(k).pos; pk=T(k).parent; ck=T(k).cost;
    if ~isnumeric(q) || ~isreal(q) || ~isvector(q) || numel(q)~=D || any(~isfinite(q(:)))
        error('plma:Position','Each node needs a finite position of the same dimension.');
    end
    if ~isnumeric(pk) || ~isreal(pk) || ~isscalar(pk) || ...
            ~isfinite(pk) || pk~=fix(pk) || pk<0 || pk>n
        error('plma:Parent','Existing parent indices must lie in 0:numel(T).');
    end
    if ~isnumeric(ck) || ~isreal(ck) || ~isscalar(ck) || ~isfinite(ck) || ck<0
        error('plma:Cost','Stored costs must be finite nonnegative scalars.');
    end
    X(k,:)=q(:)'; p(k)=pk;
end
roots=find(p==0);
if numel(roots)~=1, error('plma:Roots','The tree must have exactly one root.'); end
if p(child)==0, error('plma:RootReparent','The root cannot be reparented.'); end
if abs(T(roots).cost)>1e-10, error('plma:Cost','The root must have zero cost.'); end

state=zeros(n,1); chain=zeros(n,1);
for origin=1:n
    if state(origin)==2, continue; end
    k=origin; count=0;
    while k~=0 && state(k)==0
        state(k)=1; count=count+1; chain(count)=k; k=p(k);
    end
    if k~=0 && state(k)==1, error('plma:Cycle','Existing parent graph contains a cycle.'); end
    state(chain(1:count))=2;
end

k=parent;
while k~=0
    if k==child, error('plma:Cycle','Reparenting would create a cycle.'); end
    k=p(k);
end

T(child).parent=parent; p(child)=parent;
first_child=zeros(n,1); next_sibling=zeros(n,1);
for k=n:-1:1
    if p(k)>0
        next_sibling(k)=first_child(p(k)); first_child(p(k))=k;
    end
end
queue=zeros(n,1); queue(1)=child; head=1; tail=1; seen=false(n,1);
while head<=tail
    node=queue(head); head=head+1;
    if seen(node), error('plma:Cycle','Repeated node in the descendant graph.'); end
    seen(node)=true; pk=p(node);
    cost=T(pk).cost+norm(X(node,:)-X(pk,:));
    if ~isfinite(cost), error('plma:Cost','Recomputed geometric cost is not finite.'); end
    T(node).cost=cost;
    kid=first_child(node);
    while kid~=0
        if tail>=n, error('plma:Cycle','Descendant traversal exceeded the tree size.'); end
        tail=tail+1; queue(tail)=kid; kid=next_sibling(kid);
    end
end
end

function yes=valid_index(i,n)
yes=isnumeric(i) && isreal(i) && isscalar(i) && isfinite(i) && ...
    i==fix(i) && i>=1 && i<=n;
end
