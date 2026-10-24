function [best_d, bestApos, bestBpos, improved] = update_best_pair_by_newnode( ...
    q_new, other_nodes, best_d, bestApos, bestBpos, eps_abs)

improved = false;

if isempty(other_nodes)
    return;
end

dists = vecnorm(other_nodes - q_new, 2, 2);
[dmin, idx] = min(dists);

if dmin < best_d - eps_abs
    best_d   = dmin;
    bestApos = q_new;
    bestBpos = other_nodes(idx,:);
    improved = true;
end
end
