function r=plma_repaired(s,o)
opts=struct('max_iter',o.max_iter,'time_limit',o.time_limit);timer=tic;
[~,internal,n,it,~,~,P,tr,S]=plma_adaptive(s.start,s.goal,s.obstacles,s.space,s.sf,opts);
native_wall=toc(timer);
A=struct('X',S.A.nodes,'p',S.A.parent,'c',S.A.cost);
B=struct('X',S.B.nodes,'p',S.B.parent,'c',S.B.cost);
first=NaN;history=zeros(0,3);reason='iteration_limit';
if it<o.max_iter,reason='time_limit';end
if ~isempty(P)
 first=internal;history=[it internal sum(vecnorm(diff(P),2,2))];reason='first_connection';
end
r=struct('path',P,'raw_path',P,'processed_path',zeros(0,3),'treeA',A,'treeB',B, ...
 'iterations',it,'success',~isempty(P),'first_time',first,'elapsed',internal,'history',history);
r.details=struct('source_status','Repaired PLMA, native search parameters unchanged', ...
 'stop_reason',reason,'trace',tr,'reported_nodes',n,'costs_reconstructed',false, ...
 'options',opts,'first_path',P,'native_wall_time',native_wall);
end
