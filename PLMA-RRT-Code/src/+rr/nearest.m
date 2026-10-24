function [id,dist] = nearest(T,q)
q=q(:)';
if isempty(T.X) || numel(q)~=size(T.X,2) || any(~isfinite(q))
    error('rr:Point','Query must be finite and match a nonempty tree.');
end
d2=sum((T.X-q).^2,2); [m,id]=min(d2); dist=sqrt(m);
end
