function ok = is_valid_extension(a, b, space, obs, res, dmin)

a=a(:)'; b=b(:)';
if numel(a)~=3 || numel(b)~=3 || ~isreal([a b]) || any(~isfinite([a b]))
    ok=false; return;
end
if any(a<space(:,1)') || any(a>space(:,2)') || any(b < space(:,1)') || any(b > space(:,2)')
    ok = false;  return;
end

if collision(a, b, obs, res)
    ok = false;  return;
end

d = distanceToNearestObstacle(b, obs);
if d < dmin
    ok = false;
    return;
end
ok = true;
end
