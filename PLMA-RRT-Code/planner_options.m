function o=planner_options(id,s)
sf=s.sf;o=struct('max_iter',30000,'time_limit',120,'record_trace',false,'stop_mode','native');
switch id
 case 'GB'
  o.step=8*sf;o.goal_radius=10*sf;o.goal_bias=.1;
 case 'BI'
  o.step=Inf;o.near_radius=Inf;
 case 'F'
  o.step=Inf;o.near_radius=50*sf;o.goal_radius=10*sf;o.goal_bias=0;o.dichotomy_distance=2*sf;
 case 'Q'
  o.step=Inf;o.near_radius=50*sf;o.goal_radius=10*sf;o.goal_bias=0;o.depth=2;
 case 'BAI'
  o.step=8*sf;o.near_radius=50*sf;o.goal_radius=10*sf;
  o.safety_distance=2*sf;o.theta_max=pi;o.start_heading=[];o.goal_heading=[];
  o.postprocess_time_limit=Inf;
 case 'IBIP'
  o.step=8*sf;o.near_radius=50*sf;o.goal_radius=10*sf;o.goal_bias=.5;
  obs=s.obstacles;if isnumeric(obs),obs=num2cell(obs,2);end
  rin=cellfun(@insphere,obs);
  if isempty(rin),largest=0;else,largest=max(rin);end
  o.influence_radius=max(10*sf,largest+sqrt(3)*o.step);
  o.attraction_gain=1;o.repulsion_gain=(max(diff(s.space,1,2))/20)^2;o.goal_distance_power=2;
 case 'BIAPF'
  o.step_max=4*sf;o.near_radius=2*sf;o.connect_radius=Inf;
 case 'PLMA'
 otherwise,error('formal:Algorithm','Unknown algorithm ID.');
end
end
function r=insphere(v)
if numel(v)==4,r=v(4);elseif numel(v)==5,r=min(v(4),v(5)/2);else,r=min(v(4:6))/2;end
end
