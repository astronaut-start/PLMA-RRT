function [pathLength, elapsedTime, totalNodes, iter, avg_dist, avg_angle, path, trace, trees] = plma_adaptive(start, goal, obstacles, space, sf, opts)
if nargin<6, opts=struct(); end
if ~isfield(opts,'max_iter'), opts.max_iter=2500; end
if ~isfield(opts,'time_limit'), opts.time_limit=10; end
validateattributes(start,{'numeric'},{'vector','numel',3,'finite','real'});
validateattributes(goal,{'numeric'},{'vector','numel',3,'finite','real'});
validateattributes(space,{'numeric'},{'size',[3 2],'finite','real'});
validateattributes(sf,{'numeric'},{'scalar','positive','finite','real'});
validateattributes(opts.max_iter,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(opts.time_limit,{'numeric'},{'scalar','nonnegative','real','nonnan'});
start=start(:)'; goal=goal(:)';
if any(space(:,1)>=space(:,2)), error('plma:Workspace','Workspace bounds must have positive extent.'); end
if ~iscell(obstacles), error('plma:Obstacles','Obstacles must be a cell array.'); end
obstacles=cellfun(@(o)o(:)',obstacles,'UniformOutput',false);
input_scene=struct('space',space,'obstacles',{obstacles});
if rr.collision(start,start,input_scene)||rr.collision(goal,goal,input_scene)
    error('plma:Endpoint','Start and goal must be inside the workspace and outside closed obstacles.');
end
path=[];
avg_angle=NaN; avg_dist=NaN;


TA = init_tree(start);
TB = init_tree(goal);
iter=0;
if isequal(start,goal)
    path=start; pathLength=0; elapsedTime=0; totalNodes=2; trace=zeros(0,15);
    if nargout>=9, trees=struct('A',capture_tree_snapshot(TA),'B',capture_tree_snapshot(TB)); end
    return;
end

goal1 = goal;
goal2 = start;

p.max_iter = opts.max_iter;
p.step_size  = 8 * sf;
p.gamma      = 50 * sf;
p.conn_tol   = 10 * sf;
p.check_res  = 0.5 * sf;
p.swap_prob  = 0.5;

p.alpha_goal = 1.0;
p.alpha_rand = 0.6;
p.k_rep      = 1;
p.rep_range  = 10 * sf;
p.k_tan      = 2.0;
p.rep_eps    = 1e-3;
p.rep_clip   = 1e4;
p.use_tan    = true;

p.dmin       = 0 * sf;

p.lte_angles      = deg2rad([10 20 30 45 60 75 90]);
p.lte_pick_best   = true;
p.lte_obs_range   = max(2*p.step_size, p.rep_range);
p.fb_cone_dirs = 15;
p.fb_use_golden_phi = true;
p.fb_phi_jitter = true;

p.global_prob0      = 0.15;
p.global_prob_max   = 1;
p.global_prob_boost = 0.06;
p.global_prob_decay = 0.97;


p.stagnation_iter   = 20;
p.expand_rate       = 1.25;
p.expand_decay      = 0.96;
p.expand_max        = 3.20;

p.capsule_len_gain  = 1.10;
p.capsule_rad_gain  = 0.45;
p.capsule_len_min   = 4*p.step_size;
p.capsule_rad_min   = 2*p.step_size;
p.sample_tries      = 40;
p.progress_eps      = 1e-6;

D0 = norm(goal - start);
best_d   = D0;
bestApos = start;
bestBpos = goal;

stagn_cnt  = 0;
expand_fac = 1.0;
p_global   = p.global_prob0;

halton_k = 1;
diag_len = norm([space(1,2)-space(1,1), space(2,2)-space(2,1), space(3,2)-space(3,1)]);
p.capsule_len_max = min(diag_len, 1.5*D0);
p.capsule_rad_max = 0.50*diag_len;


p.gb_win          = 8;
p.pg_base         = 0.5;
p.pg_min          = 0;
p.pg_max          = 0.5;
p.pg_drop_k       = 10;
p.pg_recover_rate = 0.02;

pg_A = p.pg_base;  pg_B = p.pg_base;
gbA  = init_goalbias_state(p.gb_win);
gbB  = init_goalbias_state(p.gb_win);

trace=nan(p.max_iter,15); trace_count=0;
step_state=as_init_state(D0,p.step_size);
tStart = tic;


for iter = 1:p.max_iter
    if toc(tStart)>=opts.time_limit, iter=trace_count; break; end
    [q_rand, halton_k, sampInfo] = smart_sample_capsule_halton( ...
        space, start, goal, bestApos, bestBpos, best_d, expand_fac, p_global, halton_k, p);


    progress_before=step_state.progress;
    if rand < 0.5
        selected_tree=1;
        [TA,new_idx,q_new,~, pg_A, gbA, ex] = extend_star_apf(TA, space, goal1, obstacles, p, q_rand, pg_A, gbA , step_state);
        if new_idx == -1 || isempty(q_new)
        trace_count=trace_count+1;
        trace(trace_count,:)=[iter selected_tree new_idx>0 best_d progress_before ex.goal_step ex.fallback_step ex.accepted_step ex.goal_clearance ex.fallback_clearance ex.branch ex.goal_free ex.fallback_free ex.goal_open ex.fallback_open];
        step_state=as_update_state(step_state,best_d);
        continue;
        end

        nodesB = vertcat(TB.pos);
        [best_d, bestApos, bestBpos, improved] = update_best_pair_by_newnode( ...
            q_new, nodesB, best_d, bestApos, bestBpos, max(p.progress_eps*D0, 1e-9));

        TA = rewire_local(TA, new_idx, p, obstacles);

        par_idx = TA(new_idx).parent; par_pos = TA(par_idx).pos;

        if rand < 1, goal2 = q_new; end

        [found, path] = try_connect(q_new, TA, new_idx, TB, obstacles, p);

    else
        selected_tree=2;
        [TB,new_idx,q_new,~, pg_B, gbB, ex] = extend_star_apf(TB, space, goal2, obstacles, p, q_rand, pg_B, gbB , step_state);
        if new_idx == -1 || isempty(q_new)
        trace_count=trace_count+1;
        trace(trace_count,:)=[iter selected_tree new_idx>0 best_d progress_before ex.goal_step ex.fallback_step ex.accepted_step ex.goal_clearance ex.fallback_clearance ex.branch ex.goal_free ex.fallback_free ex.goal_open ex.fallback_open];
        step_state=as_update_state(step_state,best_d);
        continue;
        end

        nodesA = vertcat(TA.pos);
        [best_d, bestApos, bestBpos, improved] = update_best_pair_by_newnode( ...
            q_new, nodesA, best_d, bestApos, bestBpos, max(p.progress_eps*D0, 1e-9));

        TB = rewire_local(TB, new_idx, p, obstacles);

        par_idx = TB(new_idx).parent; par_pos = TB(par_idx).pos;

        [found, path] = try_connect(q_new, TB, new_idx, TA, obstacles, p);
    end

    trace_count=trace_count+1;
        trace(trace_count,:)=[iter selected_tree new_idx>0 best_d progress_before ex.goal_step ex.fallback_step ex.accepted_step ex.goal_clearance ex.fallback_clearance ex.branch ex.goal_free ex.fallback_free ex.goal_open ex.fallback_open];
        step_state=as_update_state(step_state,best_d);

    if improved
        stagn_cnt  = 0;
        expand_fac = max(1.0, expand_fac * p.expand_decay);
        p_global   = max(p.global_prob0, p_global * p.global_prob_decay);
    else
        stagn_cnt = stagn_cnt + 1;
        if stagn_cnt >= p.stagnation_iter
            expand_fac = min(p.expand_max, expand_fac * p.expand_rate);
            p_global   = min(p.global_prob_max, p_global + p.global_prob_boost);
            stagn_cnt  = 0;
        end
    end

    if found
        if selected_tree==2, path=flipud(path); end
        [valid_path,validation_reason]=plma_validate_path(path,start,goal,space,obstacles);
        if ~valid_path
            error('plma:InvalidSolution','Internal connection failed full-path validation: %s',validation_reason);
        end
        elapsedTime = toc(tStart);
        totalNodes  = numel(TA) + numel(TB);
        pathLength  = sum(vecnorm(diff(path),2,2));
        trace=trace(1:trace_count,:);

        avg_dist=NaN;
        avg_angle=NaN;
        if nargout >= 9
            trees = struct('A', capture_tree_snapshot(TA), ...
                'B', capture_tree_snapshot(TB));
        end
        return;
    end
end

iter=trace_count;
trace=trace(1:trace_count,:);
elapsedTime = toc(tStart);
totalNodes  = numel(TA) + numel(TB);
pathLength  = NaN; avg_dist = NaN; avg_angle = NaN;
if nargout >= 9
    trees = struct('A', capture_tree_snapshot(TA), ...
        'B', capture_tree_snapshot(TB));
end
end

function [T,new_idx,q_new,q_near, pg, gb, ex] = extend_star_apf(T, sp, goal, obs, p, q_rand, pg, gb , step_state)

new_idx = -1; q_new = []; q_near = [];
ex=struct('goal_step',NaN,'fallback_step',NaN,'accepted_step',NaN,'goal_clearance',NaN,'fallback_clearance',NaN,'branch',0,'goal_free',NaN,'fallback_free',NaN,'goal_open',NaN,'fallback_open',NaN);
nodes = vertcat(T.pos);

if rand < pg
    [q_near_goal, idx] = grnear_rand(goal, nodes);

    p_move=p; p_move.step_state=step_state;
    [ok_gb,q_new_gb,local_info]=steer_adaptive_goal(q_near_goal,goal,sp,obs,p_move);
    ex.goal_step=local_info.step; ex.goal_clearance=local_info.clearance;
    ex.goal_free=local_info.free_distance; ex.goal_open=local_info.openness;

    [pg, gb] = update_goalbias_pg(pg, gb, ok_gb, p);

    if ok_gb
        q_new = q_new_gb;
        if is_valid_extension(q_near_goal, q_new, sp, obs, p.check_res, p.dmin)
            [T,new_idx] = insert_only(T, q_new, idx, nodes, p, obs);
            q_near = q_near_goal;
            ex.accepted_step=norm(q_new-q_near); ex.branch=1;
            return;
        else
            new_idx = -1; q_new = []; q_near = [];
        end
    end
end

[qn, idx] = grnear_rand(q_rand, nodes);
p_move=p; local_state=step_state;
local_state.local_q=qn; local_state.local_clearance=as_clearance(qn,obs);
dir = apf_dir(qn, q_rand, goal, obs, p);
if ~all(dir==0)
    [p_move.step_size,local_info]=as_step(qn,obs,local_state,dir,sp);
    ex.fallback_step=p_move.step_size; ex.fallback_clearance=local_info.clearance;
    ex.fallback_free=local_info.free_distance; ex.fallback_open=local_info.openness;
    q_new = qn + p_move.step_size * dir;
    if p_move.step_size>0 && is_valid_extension(qn, q_new, sp, obs, p.check_res, p.dmin)
        [T,new_idx] = insert_only(T, q_new, idx, nodes, p, obs);
        q_near = qn;
        ex.accepted_step=norm(q_new-q_near); ex.branch=2;
        return;
    end
end

dir = q_rand - qn;
dist = norm(dir);
if dist==0, return; end

[p_move.step_size,local_info]=as_step(qn,obs,local_state,dir/dist,sp);
ex.fallback_step=p_move.step_size; ex.fallback_clearance=local_info.clearance;
    ex.fallback_free=local_info.free_distance; ex.fallback_open=local_info.openness;
q_new = qn + (p_move.step_size/dist) * dir;
if p_move.step_size>0 && is_valid_extension(qn, q_new, sp, obs, p.check_res, p.dmin)
    [T,new_idx] = insert_only(T, q_new, idx, nodes, p, obs);
    q_near = qn;
    ex.accepted_step=norm(q_new-q_near); ex.branch=3;
    return;
end

new_idx = -1; q_new = [];
end


function [T,new_idx] = insert_only(T, q_new, nn_idx, nodes, p, obs)
nNodes = numel(T);
dists  = vecnorm(nodes - q_new, 2, 2);
near   = find(dists <= p.gamma);

bestPar = nn_idx;
minCost = T(nn_idx).cost + norm(q_new - nodes(nn_idx,:));
for k = near'
    c_new = T(k).cost + norm(T(k).pos - q_new);
    if c_new < minCost && ~collision(T(k).pos, q_new, obs, p.check_res)
        bestPar = k;  minCost = c_new;
    end
end

new_idx           = nNodes + 1;
T(new_idx).pos    = q_new;
T(new_idx).parent = bestPar;
T(new_idx).cost   = minCost;
end

function T = rewire_local(T, new_idx, p, obs)
q_new  = T(new_idx).pos;
nodes  = vertcat(T.pos);
dists  = vecnorm(nodes - q_new, 2, 2);
near   = find(dists <= p.gamma);
for k = near'
    if k==new_idx || T(k).parent==0 || k==T(new_idx).parent, continue; end
    c_re = T(new_idx).cost + norm(T(k).pos - q_new);
    if c_re + eps < T(k).cost && ~collision(T(k).pos, q_new, obs, p.check_res)
        T = plma_reparent(T,k,new_idx);
    end
end
end

function [qn,idx] = grnear_rand(qr, nodes)
[~,idx] = min(vecnorm(nodes-qr,2,2));
qn = nodes(idx,:);
end

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

function [flag, path] = try_connect(q_new, TA, idxA, TB, obs, p)
flag = false; path = [];
nodesB = vertcat(TB.pos);
[~, nnB] = min(vecnorm(nodesB - q_new, 2, 2));
q_other  = nodesB(nnB,:);
if norm(q_new - q_other) <= p.conn_tol && ...
        ~collision(q_new, q_other, obs, p.check_res)
    pathA = trace_back(TA, idxA);
    pathB = trace_back(TB, nnB);
    path  = [pathA; flipud(pathB)];
    flag  = true;
end
end

function T = init_tree(root)
T(1).pos = root;
T(1).parent = 0;
T(1).cost = 0;
end


function path = trace_back(nodes, idx)
path = [];
visited=false(1,numel(nodes));
while idx ~= 0
    if ~isscalar(idx)||~isfinite(idx)||idx~=fix(idx)||idx<1||idx>numel(nodes)||visited(idx)
        error('plma:ParentChain','Invalid or cyclic parent chain while extracting path.');
    end
    visited(idx)=true;
    path = [nodes(idx).pos; path]; %#ok<AGROW>
    idx  = nodes(idx).parent;
end
end

function snapshot = capture_tree_snapshot(T)
snapshot = struct('nodes', zeros(0,3), 'parent', zeros(0,1), ...
    'edges', zeros(0,2), 'cost', zeros(0,1));
if isempty(T), return; end
snapshot.nodes = vertcat(T.pos);
snapshot.parent = reshape([T.parent], [], 1);
snapshot.cost = reshape([T.cost], [], 1);
child = find(snapshot.parent > 0);
snapshot.edges = [snapshot.parent(child), child];
end
