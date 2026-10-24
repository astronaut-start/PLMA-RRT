function [q, halton_k, info] = smart_sample_capsule_halton( ...
    space, start, goal, bestApos, bestBpos, best_d, expand_fac, p_global, halton_k, p)

info = struct('isGlobal', false, ...
              'c', [], ...
              'u', [], ...
              'halfLen', [], ...
              'r', []);

%% ===================== 1) 全局采样：Halton =====================
if rand < p_global
    q = halton_sample3(halton_k, space);
    halton_k = halton_k + 1;

    info.isGlobal = true;
    return;
end

%% ===================== 2) 局部采样：中心跳变胶囊体 =====================
if rand < 0.3
    c = 0.5 * (bestApos + bestBpos);
    dir = bestBpos - bestApos;

    if norm(dir) < 1e-12
        dir = goal - start;
    end
else
    c = 0.5 * (start + goal);
    dir = goal - start;
end

u = dir / max(norm(dir), 1e-12);

%% ===================== 3) 胶囊尺寸：随 best_d 自适应，并乘 expand_fac =====================
dref = max(best_d, p.step_size);

L = min(p.capsule_len_max, max(p.capsule_len_min, p.capsule_len_gain * dref)) * expand_fac;
r = min(p.capsule_rad_max, max(p.capsule_rad_min, p.capsule_rad_gain * dref)) * expand_fac;
halfLen = 0.5 * L;

%% ===================== 4) 在胶囊体体积内采样，并确保在边界内 =====================
q = sample_point_in_capsule(space, c, u, halfLen, r, p.sample_tries);

if ~in_space(q, space)
    q = halton_sample3(halton_k, space);
    halton_k = halton_k + 1;
end

info.isGlobal = false;
info.c = c;
info.u = u;
info.halfLen = halfLen;
info.r = r;
end

function tf = in_space(q, space)
tf = (q(1) >= space(1,1) && q(1) <= space(1,2) && ...
      q(2) >= space(2,1) && q(2) <= space(2,2) && ...
      q(3) >= space(3,1) && q(3) <= space(3,2));
end
