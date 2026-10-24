function free = as_free_distance(q, d, obstacles, space, max_range)

validateattributes(q, {'numeric'}, {'real', 'finite', 'vector', 'numel', 3}, ...
    mfilename, 'q');
validateattributes(d, {'numeric'}, {'real', 'finite', 'vector', 'numel', 3}, ...
    mfilename, 'd');
validateattributes(space, {'numeric'}, {'real', 'finite', 'size', [3 2]}, ...
    mfilename, 'space');
validateattributes(max_range, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'nonnegative'}, mfilename, 'max_range');
q = double(q(:).');
d = double(d(:).');
direction_length = norm(d);
if direction_length == 0
    error('as_free_distance:ZeroDirection', '射线方向不得为零。');
end
d = d/direction_length;
if any(space(:,1) > space(:,2))
    error('as_free_distance:InvalidSpace', '空间各维下界不得大于上界。');
end
if ~isempty(obstacles) && ~iscell(obstacles)
    error('as_free_distance:ObstacleContainer', '障碍集合必须是 cell 数组。');
end
free = double(max_range);
if any(q < space(:,1).') || any(q > space(:,2).')
    free = 0;
    return;
end

for axis = 1:3
    if d(axis) > 0
        exit_distance = (space(axis,2)-q(axis))/d(axis);
    elseif d(axis) < 0
        exit_distance = (space(axis,1)-q(axis))/d(axis);
    else
        continue;
    end
    free = min(free, max(0, exit_distance));
end
if free == 0
    return;
end

for k = 1:numel(obstacles)
    o = obstacles{k};
    if ~isnumeric(o) || ~isreal(o) || ~isvector(o) || ...
            any(~isfinite(o(:))) || ~ismember(numel(o), [4 5 6])
        error('as_free_distance:ObstacleFormat', ...
            '第 %d 个障碍必须是有限实数的 4、5 或 6 元向量。', k);
    end
    o = double(o(:).');
    if any(o(4:end) < 0)
        error('as_free_distance:ObstacleSize', '第 %d 个障碍的尺寸不得为负。', k);
    end
    x = q-o(1:3);
    switch numel(o)
        case 4
            c = dot(x,x)-o(4)^2;
            if c <= 0
                free = 0;
                return;
            end
            [entry, leave] = radial_interval(dot(d,d), dot(x,d), c);
        case 6
            half = o(4:6)/2;
            if all(abs(x) <= half)
                free = 0;
                return;
            end
            entry = -inf;
            leave = inf;
            for axis = 1:3
                [lo, hi] = slab_interval(x(axis), d(axis), half(axis));
                entry = max(entry,lo);
                leave = min(leave,hi);
            end
        case 5
            radial_c = dot(x(1:2),x(1:2))-o(4)^2;
            if radial_c <= 0 && abs(x(3)) <= o(5)/2
                free = 0;
                return;
            end
            [radial_entry, radial_leave] = radial_interval( ...
                dot(d(1:2),d(1:2)), dot(x(1:2),d(1:2)), radial_c);
            [axial_entry, axial_leave] = slab_interval(x(3),d(3),o(5)/2);
            entry = max(radial_entry,axial_entry);
            leave = min(radial_leave,axial_leave);
    end
    first = max(0,entry);
    if first <= min(free,leave)
        free = first;
    end
end
end

function [entry, leave] = slab_interval(x, d, half)
if d == 0
    if abs(x) <= half
        entry = -inf;
        leave = inf;
    else
        entry = inf;
        leave = -inf;
    end
else
    ends = [(-half-x)/d, (half-x)/d];
    entry = min(ends);
    leave = max(ends);
end
end

function [entry, leave] = radial_interval(a,b,c)
if a == 0
    if c <= 0
        entry = -inf;
        leave = inf;
    else
        entry = inf;
        leave = -inf;
    end
    return;
end
discriminant = b*b-a*c;
roundoff = 32*eps(max(abs([b*b,a*c])));
if discriminant < -roundoff
    entry = inf;
    leave = -inf;
    return;
end
root = sqrt(max(0,discriminant));
if b >= 0
    numerator = -b-root;
else
    numerator = -b+root;
end
if numerator == 0
    entry = -b/a;
    leave = entry;
else
    roots = [numerator/a, c/numerator];
    entry = min(roots);
    leave = max(roots);
end
end
