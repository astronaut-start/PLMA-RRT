function [T,id] = insert(T,q,parent)
q=q(:)'; n=size(T.X,1);
if numel(q)~=size(T.X,2) || any(~isfinite(q))
    error('rr:Point','Inserted point must be finite and match the tree dimension.');
end
if ~isscalar(parent) || ~isfinite(parent) || parent~=fix(parent) || parent<1 || parent>n
    error('rr:Parent','Parent must index an existing node.');
end
id=n+1; T.X(id,:)=q; T.p(id,1)=parent;
T.c(id,1)=T.c(parent)+norm(q-T.X(parent,:));
end
