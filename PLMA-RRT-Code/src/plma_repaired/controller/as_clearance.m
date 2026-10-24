function clearance = as_clearance(q, obstacles)

validateattributes(q, {'numeric'}, {'real', 'finite', 'vector', 'numel', 3}, ...
    mfilename, 'q');
q = double(q(:).');
clearance = inf;
if isempty(obstacles)
    return;
end
if ~iscell(obstacles)
    error('as_clearance:ObstacleContainer', '障碍集合必须是 cell 数组。');
end

for k = 1:numel(obstacles)
    o = obstacles{k};
    if ~isnumeric(o) || ~isreal(o) || ~isvector(o) || ...
            any(~isfinite(o(:))) || ~ismember(numel(o), [4 5 6])
        error('as_clearance:ObstacleFormat', ...
            '第 %d 个障碍必须是有限实数的 4、5 或 6 元向量。', k);
    end
    o = double(o(:).');
    if any(o(4:end) < 0)
        error('as_clearance:ObstacleSize', '第 %d 个障碍的尺寸不得为负。', k);
    end
    relative = q - o(1:3);
    switch numel(o)
        case 4
            distance = max(norm(relative) - o(4), 0);
        case 6
            outside = max(abs(relative) - o(4:6)/2, 0);
            distance = norm(outside);
        case 5
            radial = max(hypot(relative(1), relative(2)) - o(4), 0);
            axial = max(abs(relative(3)) - o(5)/2, 0);
            distance = hypot(radial, axial);
    end
    clearance = min(clearance, distance);
end
end
