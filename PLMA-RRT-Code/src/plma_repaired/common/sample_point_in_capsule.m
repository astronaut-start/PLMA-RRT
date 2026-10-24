function q = sample_point_in_capsule(space, c, u, halfLen, r, maxTries)

[v, w] = orthonormal_basis(u);

for k = 1:maxTries
    s = (2*rand - 1) * (halfLen + r);
    a = abs(s);

    if a <= halfLen
        rhoMax = r;
    else
        dcap = a - halfLen;
        if dcap > r
            continue;
        end
        rhoMax = sqrt(max(0, r^2 - dcap^2));
    end

    rho   = rhoMax * sqrt(rand);
    theta = 2*pi*rand;

    offset = rho*cos(theta)*v + rho*sin(theta)*w;

    q = c + s*u + offset;

    if in_space(q, space)
        return;
    end
end

q = clamp_to_space(c, space);
end

function tf = in_space(q, space)
tf = (q(1) >= space(1,1) && q(1) <= space(1,2) && ...
      q(2) >= space(2,1) && q(2) <= space(2,2) && ...
      q(3) >= space(3,1) && q(3) <= space(3,2));
end

function q = clamp_to_space(q, space)
q(1) = min(max(q(1), space(1,1)), space(1,2));
q(2) = min(max(q(2), space(2,1)), space(2,2));
q(3) = min(max(q(3), space(3,1)), space(3,2));
end
