function T = reparent(T,child,parent)
n=size(T.X,1);
if ~isscalar(child) || ~isscalar(parent) || ...
        any(~isfinite([child parent])) || any([child parent]~=fix([child parent])) || ...
        any([child parent]<1) || any([child parent]>n)
    error('rr:Parent','Child and parent must index existing nodes.');
end
if T.p(child)==0, error('rr:RootReparent','The root cannot be reparented.'); end
visited=false(n,1); k=parent;
while k~=0
    if ~isfinite(k) || k<1 || k>n || k~=fix(k), error('rr:Parent','Malformed parent chain.'); end
    if k==child || visited(k), error('rr:Cycle','Reparenting would create or encounter a cycle.'); end
    visited(k)=true; k=T.p(k);
end
T.p(child)=parent;
queue=zeros(n,1); queue(1)=child; head=1; tail=1; visited=false(n,1);
while head<=tail
    k=queue(head); head=head+1;
    if visited(k), error('rr:Cycle','Malformed descendant graph.'); end
    visited(k)=true; pk=T.p(k);
    T.c(k)=T.c(pk)+norm(T.X(k,:)-T.X(pk,:));
    kids=find(T.p==k);
    if tail+numel(kids)>n, error('rr:Cycle','Malformed descendant graph.'); end
    queue(tail+(1:numel(kids)))=kids; tail=tail+numel(kids);
end
end
