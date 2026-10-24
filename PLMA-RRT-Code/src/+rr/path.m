function P = path(T,id)
n=size(T.X,1);
if ~isscalar(id) || ~isfinite(id) || id~=fix(id) || id<1 || id>n
    error('rr:Node','Target node must be a valid index.');
end
indices=zeros(n,1); visited=false(n,1); count=0;
while id~=0
    if ~isfinite(id) || id<1 || id>n || id~=fix(id), error('rr:Parent','Malformed parent index.'); end
    if visited(id), error('rr:Cycle','Tree contains a parent cycle.'); end
    visited(id)=true; count=count+1; indices(count)=id; id=T.p(id);
end
P=T.X(flipud(indices(1:count)),:);
end
