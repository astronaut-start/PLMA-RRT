function [dist, closest_point] = distanceToNearestObstacle(point, obstacles)

dist = inf;
closest_point = [NaN, NaN, NaN];

px = point(1);
py = point(2);
pz = point(3);

for i = 1:numel(obstacles)
    obs = obstacles{i};

    if length(obs) == 6
        cx = obs(1); cy = obs(2); cz = obs(3);
        dx = obs(4); dy = obs(5); dz = obs(6);

        x_min = cx - dx/2; x_max = cx + dx/2;
        y_min = cy - dy/2; y_max = cy + dy/2;
        z_min = cz - dz/2; z_max = cz + dz/2;

        closest_x = clamp(px, x_min, x_max);
        closest_y = clamp(py, y_min, y_max);
        closest_z = clamp(pz, z_min, z_max);
        curPoint = [closest_x, closest_y, closest_z];

        dxp = abs(px - cx) - dx/2;
        dyp = abs(py - cy) - dy/2;
        dzp = abs(pz - cz) - dz/2;
        dxp = max(dxp, 0); dyp = max(dyp, 0); dzp = max(dzp, 0);
        curDist = sqrt(dxp^2 + dyp^2 + dzp^2);

    elseif length(obs) == 4
        cx = obs(1); cy = obs(2); cz = obs(3);
        r = obs(4);

        vec = [px - cx, py - cy, pz - cz];
        center_dist = norm(vec);

        if center_dist == 0
            curPoint = [cx + r, cy, cz];
        else
            scale = r / center_dist;
            curPoint = [cx, cy, cz] + vec * scale;
        end

        curDist = max(center_dist - r, 0);

    elseif length(obs) == 5
        cx = obs(1); cy = obs(2); cz = obs(3);
        r  = obs(4); h  = obs(5);

        z_min = cz - h/2;
        z_max = cz + h/2;

        dx = px - cx;
        dy = py - cy;
        dist_xy = sqrt(dx^2 + dy^2);

        if dist_xy > 0
            scale = r / dist_xy;
            closest_x = cx + dx * scale;
            closest_y = cy + dy * scale;
        else
            closest_x = cx + r;
            closest_y = cy;
        end

        if pz < z_min
            closest_z = z_min;
        elseif pz > z_max
            closest_z = z_max;
        else
            closest_z = pz;
        end

        vertical_dist = 0;
        if pz < z_min
            vertical_dist = z_min - pz;
        elseif pz > z_max
            vertical_dist = pz - z_max;
        end

        if vertical_dist > 0
            dist_side = norm([px - closest_x, py - closest_y, pz - closest_z]);

            if dist_xy <= r
                dist_top_bottom = vertical_dist;
            else
                dist_top_bottom = sqrt((dist_xy - r)^2 + vertical_dist^2);
            end

            if dist_top_bottom < dist_side
                closest_z = z_min * (pz < z_min) + z_max * (pz > z_max);
                if dist_xy > 0
                    closest_x = cx + dx * r/dist_xy;
                    closest_y = cy + dy * r/dist_xy;
                end
            end
        end
        curPoint = [closest_x, closest_y, closest_z];
        curDist = norm(point - curPoint);

    else
        error('未知障碍物格式 (长度=%d)', length(obs));
    end

    if curDist < dist
        dist = curDist;
        closest_point = curPoint;
    end
end
end

function v = clamp(value, min_val, max_val)
v = max(min(value, max_val), min_val);
end
