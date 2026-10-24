function T = new_tree(root)
root=root(:)';
if ~ismember(numel(root),[2 3]) || any(~isfinite(root))
    error('rr:Root','Root must be a finite point in 2D or 3D.');
end
T=struct('X',root,'p',0,'c',0);
end
