function [ok,q_best,selected_info] = steer_adaptive_goal(q_from,goal,sp,obs,p)

ok = false;
q_best = [];
local_state=p.step_state;
local_state.local_q=q_from;local_state.local_clearance=as_clearance(q_from,obs);
selected_info=struct('step',0,'clearance',local_state.local_clearance,'free_distance',NaN,'openness',NaN);

dir  = goal - q_from;
dist = norm(dir);
if dist < 1e-12
    return;
end

d0    = dir / dist;
[h,step_info]=as_step(q_from,obs,local_state,d0,sp);
selected_info=step_info;
q_try = q_from + h * d0;

if h>0 && is_valid_extension(q_from, q_try, sp, obs, p.check_res, p.dmin)
    ok = true;
    q_best = q_try;
    return;
end

bestScore = inf;

[n, dminObs] = nearest_obstacle_normal(q_from, obs);

useTrueLTE = (norm(n) > 1e-12) && (dminObs <= p.lte_obs_range);

if useTrueLTE
    n = normalize_vec(n);

    t1 = cross(n, d0);

    if norm(t1) < 1e-9
        [vtmp, ~] = orthonormal_basis(d0);
        t1 = cross(n, vtmp);
    end
    t1 = normalize_vec(t1);

    t2 = normalize_vec(cross(n, t1));

    for a = p.lte_angles
        ca = cos(a);
        sa = sin(a);

        cand_dirs = [
            ca*d0 + sa*t1;
            ca*d0 - sa*t1;
            ca*d0 + sa*t2;
            ca*d0 - sa*t2;
            ca*d0 + sa*normalize_vec(t1+t2);
            ca*d0 + sa*normalize_vec(t1-t2);
        ];

        for k = 1:size(cand_dirs,1)
            d = normalize_vec(cand_dirs(k,:));
            if norm(d)<1e-12,continue;end
            [h,step_info]=as_step(q_from,obs,local_state,d,sp);
            q = q_from + h * d;

            if h>0 && is_valid_extension(q_from, q, sp, obs, p.check_res, p.dmin)
                if ~p.lte_pick_best
                    ok = true; q_best = q; selected_info=step_info; return;
                else
                    score = norm(goal - q);
                    if score < bestScore
                        bestScore = score;
                        q_best = q; selected_info=step_info;
                        ok = true;
                    end
                end
            end
        end

    end

    if ok, return; end
end

[e1, e2] = cone_frame_from_axis(d0);

phi0 = 0;
if p.fb_phi_jitter
    phi0 = 2*pi*rand;
end

for a = p.lte_angles
    ca = cos(a);
    sa = sin(a);

    for m = 1:p.fb_cone_dirs
        if p.fb_use_golden_phi
            golden = 2.399963229728653;
            phi = mod((m-1)*golden + phi0, 2*pi);
        else
            phi = (m-1) * (2*pi/p.fb_cone_dirs) + phi0;
        end

        t = cos(phi)*e1 + sin(phi)*e2;

        d = normalize_vec(ca*d0 + sa*t);
        if norm(d)<1e-12,continue;end

        [h,step_info]=as_step(q_from,obs,local_state,d,sp);
            q = q_from + h * d;

        if h>0 && is_valid_extension(q_from, q, sp, obs, p.check_res, p.dmin)
            if ~p.lte_pick_best
                ok = true; q_best = q; selected_info=step_info; return;
            else
                score = norm(goal - q);
                if score < bestScore
                    bestScore = score;
                    q_best = q; selected_info=step_info;
                    ok = true;
                end
            end
        end
    end
end

end

function [n_best, d_best] = nearest_obstacle_normal(q, obstacles)

n_best = [0 0 0];
d_best = inf;

for i = 1:numel(obstacles)
    o = obstacles{i};
    if ~isnumeric(o) || isempty(o), continue; end

    if numel(o) == 4
        [n, d] = closest_sphere(q, o);

    elseif numel(o) == 6
        [n, d] = closest_box(q, o);

    elseif numel(o) == 5
        [n, d] = closest_cylinder_z(q, o);

    else
        continue;
    end

    if d < d_best
        d_best = d;
        n_best = n;
    end
end
end

function [n, d] = closest_sphere(q, o)
c = o(1:3);
r = o(4);

v  = q - c;
nv = norm(v);

if nv < 1e-12
    n  = [1 0 0];
    cp = c + r*n;
else
    n  = v / nv;
    cp = c + r*n;
end

d = norm(q - cp);
end

function [n, d] = closest_box(q, o)
c = o(1:3);
h = o(4:6)/2;
bmin = c - h;
bmax = c + h;

cp = min(max(q, bmin), bmax);

inside = all(q >= bmin) && all(q <= bmax);

if ~inside
    v  = q - cp;
    nv = norm(v);
    if nv < 1e-12
        n = [1 0 0];
    else
        n = v / nv;
    end
else
    dx = h(1) - abs(q(1)-c(1));
    dy = h(2) - abs(q(2)-c(2));
    dz = h(3) - abs(q(3)-c(3));
    [~, id] = min([dx dy dz]);

    cp = q;
    if id == 1
        s = sign(q(1)-c(1)); if s==0, s=1; end
        n = [s 0 0];
        cp(1) = c(1) + s*h(1);
    elseif id == 2
        s = sign(q(2)-c(2)); if s==0, s=1; end
        n = [0 s 0];
        cp(2) = c(2) + s*h(2);
    else
        s = sign(q(3)-c(3)); if s==0, s=1; end
        n = [0 0 s];
        cp(3) = c(3) + s*h(3);
    end
end

d = norm(q - cp);
end

function [n, d] = closest_cylinder_z(q, o)
cx=o(1); cy=o(2); cz=o(3); r=o(4); h=o(5);
zmin = cz - h/2;
zmax = cz + h/2;

dx  = q(1) - cx;
dy  = q(2) - cy;
rho = sqrt(dx*dx + dy*dy);

if rho < 1e-12
    ur = [1 0 0];
else
    ur = [dx dy 0]/rho;
end

insideRad = (rho <= r);
insideZ   = (q(3) >= zmin) && (q(3) <= zmax);
inside    = insideRad && insideZ;

if ~inside
    zc = min(max(q(3), zmin), zmax);

    if insideZ && (rho > r)
        cp = [cx cy 0] + r*ur + [0 0 zc];
        n  = ur;

    elseif (~insideZ) && insideRad
        zcap = zmin; if q(3) > zmax, zcap = zmax; end
        cp = [q(1) q(2) zcap];
        n  = [0 0 sign(q(3)-cz)];
        if n(3)==0, n=[0 0 1]; end

    else
        zcap = zmin; if q(3) > zmax, zcap = zmax; end
        cp = [cx cy 0] + r*ur + [0 0 zcap];
        v = q - cp;
        n = normalize_vec(v);
        if norm(n) < 1e-12, n = [0 0 1]; end
    end
else
    ds = r - rho;
    dc = h/2 - abs(q(3)-cz);

    if ds <= dc
        cp = [cx cy q(3)] + r*ur;
        n  = ur;
    else
        if q(3) >= cz
            cp = [q(1) q(2) zmax];
            n  = [0 0 1];
        else
            cp = [q(1) q(2) zmin];
            n  = [0 0 -1];
        end
    end
end

d = norm(q - cp);
end

function x = normalize_vec(x)
n = norm(x);
if n < 1e-12
    x = [0 0 0];
else
    x = x / n;
end
end

function [v,w] = orthonormal_basis(u)
u = u(:)' / max(norm(u),1e-12);

if abs(u(1)) < 0.9
    a = [1 0 0];
else
    a = [0 1 0];
end

v = cross(u,a); v = v / max(norm(v),1e-12);
w = cross(u,v); w = w / max(norm(w),1e-12);
end


function [e1, e2] = cone_frame_from_axis(d0)
d0 = normalize_vec(d0);

a = [1 0 0];
if abs(dot(a, d0)) > 0.9
    a = [0 1 0];
end

e1 = a - dot(a, d0)*d0;
e1 = normalize_vec(e1);

e2 = cross(d0, e1);
e2 = normalize_vec(e2);

if norm(e1) < 1e-12 || norm(e2) < 1e-12
    a = [0 0 1];
    if abs(dot(a, d0)) > 0.9
        a = [0 1 0];
    end
    e1 = normalize_vec(a - dot(a,d0)*d0);
    e2 = normalize_vec(cross(d0, e1));
end
end
