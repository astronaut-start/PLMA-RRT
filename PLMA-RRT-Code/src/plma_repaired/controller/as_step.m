function [h, info] = as_step(q, obstacles, s, dir, space)

narginchk(5,5);
if isfield(s,'local_q') && isfield(s,'local_clearance') && ...
        isequal(q(:).',s.local_q)
    clearance = s.local_clearance;
else
    clearance = as_clearance(q, obstacles);
end
clearance_ratio = min(1, max(0, clearance/s.h0));
free = as_free_distance(q,dir,obstacles,space,s.hmax);
openness = 0.25*clearance_ratio + 0.75*(free/s.hmax);
progress = s.progress;
desired = s.hmin + (s.h0*(1 + 0.75*progress)-s.hmin)*openness;
desired = min(s.hmax,max(s.hmin,desired));
if free < s.hmax
    h = min(desired,0.90*free);
else
    h = desired;
end
blocked = h < s.hmin;
if blocked
    h = 0;
end
info = struct('progress', progress, 'clearance', clearance, ...
    'openness', openness, 'step', h, 'free_distance', free, ...
    'desired_step', desired, 'blocked', blocked);
end
