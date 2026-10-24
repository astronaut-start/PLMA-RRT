function s = as_update_state(s, new_best_d)

validateattributes(new_best_d, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'nonnegative'}, mfilename, 'new_best_d');

next_best = min(s.best_d, double(new_best_d));
drop = max(0, s.best_d - next_best);
s.best_d = next_best;
s.window_gains = [s.window_gains, drop];
if numel(s.window_gains) > s.W
    s.window_gains = s.window_gains(end-s.W+1:end);
end
s.attempt_count = s.attempt_count + 1;
s.last_drop = drop;
s.mean_drop = mean(s.window_gains);

s.progress = min(1, max(0, s.mean_drop / s.progress_scale));
end
