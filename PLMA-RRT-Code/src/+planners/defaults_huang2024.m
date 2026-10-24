function o = defaults_huang2024(s, overrides)
span=s.space(:,2)-s.space(:,1); L=norm(span);
if L<=0, error('huang:Scene','The scene must have positive extent.'); end
o=struct('max_iter',1000,'time_limit',Inf,'stop_mode','native', ...
    'step_max',10,'near_radius',5,'connect_radius',Inf, ...
    'Katt',1,'Krep',1,'rep_power',2,'influence_radius',0.1*L, ...
    'goal_region_radius',0.1*L,'density_alpha',0.5,'density_beta',0.5, ...
    'goal_bias_initial',0.2,'goal_bias_distance_decay',1/L, ...
    'goal_bias_iteration_decay',0.001,'step_density_decay',0.1, ...
    'target_switch_probability',0.1,'region_sigma',0.15*L, ...
    'record_trace',false);
if nargin>1
    keys=fieldnames(overrides);
    for k=1:numel(keys)
        if ~isfield(o,keys{k}), error('huang:Option','Unknown option: %s',keys{k}); end
        o.(keys{k})=overrides.(keys{k});
    end
end
if ~any(strcmp(o.stop_mode,{'native','first_solution','first_valid'}))
    error('huang:NativeFirstSolution', ...
        'Huang Figure 3 returns at first connection; use native or first_solution (legacy first_valid also accepted).');
end
validateattributes(o.max_iter,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(o.time_limit,{'numeric'},{'scalar','positive'});
validateattributes(o.connect_radius,{'numeric'},{'scalar','positive'});
positive={'step_max','near_radius','influence_radius', ...
    'goal_region_radius','region_sigma'};
for k=1:numel(positive)
    validateattributes(o.(positive{k}),{'numeric'},{'scalar','finite','positive'});
end
nonnegative={'Katt','Krep','density_alpha','goal_bias_distance_decay', ...
    'goal_bias_iteration_decay','step_density_decay'};
for k=1:numel(nonnegative)
    validateattributes(o.(nonnegative{k}),{'numeric'},{'scalar','finite','nonnegative'});
end
validateattributes(o.rep_power,{'numeric'},{'scalar','finite','>=',1});
prob={'density_beta','goal_bias_initial','target_switch_probability'};
for k=1:numel(prob)
    validateattributes(o.(prob{k}),{'numeric'},{'scalar','finite','>=',0,'<=',1});
end
validateattributes(o.record_trace,{'logical'},{'scalar'});
end
