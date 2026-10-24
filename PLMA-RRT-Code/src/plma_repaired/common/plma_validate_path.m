function [valid,reason]=plma_validate_path(P,start,goal,space,obstacles)
valid=false; reason='';
if ~isnumeric(P)||isempty(P)||size(P,2)~=3||any(~isfinite(P),'all')
    reason='Empty, nonfinite or dimensionally invalid path.'; return;
end
if norm(P(1,:)-start(:)')>1e-9 || norm(P(end,:)-goal(:)')>1e-9
    reason='Endpoints do not match start and goal in order.'; return;
end
s=struct('space',space,'obstacles',{obstacles});
for k=1:size(P,1)
    if rr.collision(P(k,:),P(k,:),s)
        reason=sprintf('Node %d is occupied or outside workspace.',k); return;
    end
    if k>1 && rr.collision(P(k-1,:),P(k,:),s)
        reason=sprintf('Complete path segment %d intersects an obstacle or leaves workspace.',k-1); return;
    end
end
valid=true;
end
