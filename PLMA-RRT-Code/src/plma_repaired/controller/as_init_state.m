function s = as_init_state(D0, h0)

validateattributes(D0, {'numeric'}, {'real', 'finite', 'scalar', 'nonnegative'}, ...
    mfilename, 'D0');
validateattributes(h0, {'numeric'}, {'real', 'finite', 'scalar', 'positive'}, ...
    mfilename, 'h0');

s.D0 = double(D0);
s.best_d = double(D0);
s.h0 = double(h0);
s.hmin = 0.35 * s.h0;
s.hmax = 1.75 * s.h0;
s.W = 12;
s.progress_scale = 0.25 * s.h0;
s.window_gains = zeros(1, 0);
s.attempt_count = 0;
s.last_drop = 0;
s.mean_drop = 0;
s.progress = 0.5;
end
