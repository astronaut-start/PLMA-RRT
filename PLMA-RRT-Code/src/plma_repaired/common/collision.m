function flag = collision(a, b, obs, res) %#ok<INUSD>
narginchk(3,4);
if ~isnumeric(a) || ~isreal(a) || ~isvector(a) || numel(a)~=3 || ...
        ~isnumeric(b) || ~isreal(b) || ~isvector(b) || numel(b)~=3
    error('plma:CollisionDimension','Endpoints must be real numeric 3D vectors.');
end
if ~iscell(obs)
    error('plma:CollisionObstacles','Obstacles must be a cell array.');
end
a=a(:)'; b=b(:)';
if any(~isfinite(a)) || any(~isfinite(b))
    flag=true; return;
end
flag=false;
for k=1:numel(obs)
    obstacle=obs{k};
    if ~isnumeric(obstacle) || ~isreal(obstacle) || ~isvector(obstacle) || ...
            ~ismember(numel(obstacle),[4 5 6])
        error('plma:CollisionObstacle','Each obstacle must be a real numeric sphere, box or cylinder vector.');
    end
    if rr.obstacle_hit(a,b,obstacle)
        flag=true; return;
    end
end
end
