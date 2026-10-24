function hit = collision(a,b,s)
a=a(:)'; b=b(:)'; D=numel(a);
if numel(b)~=D || ~ismember(D,[2 3])
    error('rr:Dimension','Endpoints must have the same dimension (2 or 3).');
end
if ~isfield(s,'space') || ~isequal(size(s.space),[D 2]) || ...
        any(~isfinite(s.space),'all') || any(s.space(:,1)>s.space(:,2))
    error('rr:Scene','scene.space must contain finite ordered bounds in D-by-2 form.');
end
if ~isfield(s,'obstacles') || ~iscell(s.obstacles)
    error('rr:Scene','scene.obstacles must be a cell array.');
end
hit=any(~isfinite(a)) || any(~isfinite(b)) || ...
    any(a<s.space(:,1)') || any(a>s.space(:,2)') || ...
    any(b<s.space(:,1)') || any(b>s.space(:,2)');
if hit, return; end
for k=1:numel(s.obstacles)
    if rr.obstacle_hit(a,b,s.obstacles{k}), hit=true; return; end
end
end
