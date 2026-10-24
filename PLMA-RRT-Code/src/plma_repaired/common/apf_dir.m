function [dir, F_attr, F_rep, F_tan] = apf_dir(q_near, q_rand, goal, obs, p)

%% ---------- 吸引力 ----------
F_attr = p.alpha_goal * unit(goal  - q_near) + ...
    p.alpha_rand * unit(q_rand - q_near);

att_dir = unit(F_attr);

%% ---------- 斥力 + 切向 ----------
F_rep = zeros(1,3);
F_tan = zeros(1,3);

for j = 1:numel(obs)
    o = obs{j};

    if numel(o) == 4
        [n, d] = closest_sphere(q_near, o);
    elseif numel(o) == 6
        [n, d] = closest_box(q_near, o);
    elseif numel(o) == 5
        [n, d] = closest_cylinder_z(q_near, o);
    else
        continue;
    end

    if d >= p.rep_range, continue; end

    d_eff = max(d, p.rep_eps);
    n = unit(n);

    s = (p.rep_range - d) / p.rep_range;
    s = max(0, min(1, s));

    mag = p.k_rep * (1/d_eff - 1/p.rep_range) / (d_eff^2);
    mag = mag * (s^2);

    mag = min(mag, p.rep_clip);

    F_rep = F_rep + mag * n;

    if p.use_tan
        if dot(att_dir, n) < 0
            t = att_dir - dot(att_dir, n) * n;
            if norm(t) > 1e-9
                t = unit(t);
                F_tan = F_tan + (p.k_tan * mag) * t;
            end
        end
    end
end

%% ---------- 合力 ----------
F = F_attr + F_rep + F_tan;

if all(abs(F) < 1e-12)
    dir = [0 0 0];
else
    dir = unit(F);
end
end

function [n, d] = closest_sphere(q, o)
c = o(1:3); r = o(4);
v = q - c; nv = norm(v);
if nv < 1e-12
    n = [1 0 0];
    cp = c + r*n;
else
    n = v / nv;
    cp = c + r*n;
end
d = norm(q - cp);
end

function [n, d] = closest_box(q, o)
c = o(1:3); h = o(4:6)/2;
bmin = c - h; bmax = c + h;

cp = min(max(q, bmin), bmax);
inside = all(q >= bmin) && all(q <= bmax);

if ~inside
    v = q - cp;
    if norm(v) < 1e-12, n = [1 0 0]; else, n = v / norm(v); end
else
    dx = h(1) - abs(q(1)-c(1));
    dy = h(2) - abs(q(2)-c(2));
    dz = h(3) - abs(q(3)-c(3));
    [~,id] = min([dx dy dz]);

    cp = q;
    if id == 1
        s = sign(q(1)-c(1)); if s==0, s=1; end
        n = [s 0 0]; cp(1) = c(1) + s*h(1);
    elseif id == 2
        s = sign(q(2)-c(2)); if s==0, s=1; end
        n = [0 s 0]; cp(2) = c(2) + s*h(2);
    else
        s = sign(q(3)-c(3)); if s==0, s=1; end
        n = [0 0 s]; cp(3) = c(3) + s*h(3);
    end
end

d = norm(q - cp);
end

function [n, d] = closest_cylinder_z(q, o)
cx=o(1); cy=o(2); cz=o(3); r=o(4); h=o(5);
zmin = cz - h/2; zmax = cz + h/2;

dx = q(1)-cx; dy = q(2)-cy;
rho = sqrt(dx*dx + dy*dy);
if rho < 1e-12, ur = [1 0 0]; else, ur = [dx dy 0]/rho; end

insideRad = (rho <= r);
insideZ   = (q(3) >= zmin) && (q(3) <= zmax);
inside = insideRad && insideZ;

if ~inside
    zc = min(max(q(3), zmin), zmax);

    if insideZ && (rho > r)
        cp = [cx cy 0] + r*ur + [0 0 zc];
        n  = ur;
    elseif (~insideZ) && insideRad
        zcap = zmin; if q(3) > zmax, zcap = zmax; end
        cp = [q(1) q(2) zcap];
        n  = [0 0 sign(q(3)-cz)]; if n(3)==0, n=[0 0 1]; end
    else
        zcap = zmin; if q(3) > zmax, zcap = zmax; end
        cp = [cx cy 0] + r*ur + [0 0 zcap];
        v = q - cp;
        if norm(v)<1e-12, n=[0 0 1]; else, n=v/norm(v); end
    end
else
    ds = r - rho;
    dc = h/2 - abs(q(3)-cz);
    if ds <= dc
        cp = [cx cy q(3)] + r*ur;
        n  = ur;
    else
        if q(3) >= cz, cp=[q(1) q(2) zmax]; n=[0 0 1];
        else,          cp=[q(1) q(2) zmin]; n=[0 0 -1];
        end
    end
end

d = norm(q - cp);
end

function u = unit(v)
n = norm(v);
if n < 1e-12, u = [0 0 0]; else, u = v / n; end
end
