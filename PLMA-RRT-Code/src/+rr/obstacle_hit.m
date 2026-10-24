function hit = obstacle_hit(a, b, o)
a = a(:)'; b = b(:)'; o = o(:)'; D = numel(a);
if ~ismember(D,[2 3]) || numel(b) ~= D
    error('rr:Dimension','Endpoints must have the same dimension (2 or 3).');
end
if any(~isfinite(a)) || any(~isfinite(b))
    hit = true; return;
end
if any(~isfinite(o))
    error('rr:Obstacle','Obstacle entries must be finite.');
end
if D == 3 && numel(o) == 4 || D == 2 && numel(o) == 3
    if o(end) < 0, error('rr:Obstacle','Radius must be nonnegative.'); end
    v = b-a; w = a-o(1:D); denom = dot(v,v);
    if denom == 0, t = 0; else, t = min(1,max(0,-dot(w,v)/denom)); end
    tol = 64*eps(max([1,abs(a),abs(b),abs(o)]));
    hit = norm(w+t*v) <= o(end) + tol;
elseif numel(o) == 2*D
    widths = o(D+1:end);
    if any(widths < 0), error('rr:Obstacle','Box widths must be nonnegative.'); end
    lo = o(1:D)-widths/2; hi = o(1:D)+widths/2;
    tol = 64*eps(max([1,abs(a),abs(b),abs(o)]));
    tlo = 0; thi = 1; v = b-a;
    for j = 1:D
        if v(j) == 0
            if a(j) < lo(j)-tol || a(j) > hi(j)+tol, hit=false; return; end
        else
            ts = sort([(lo(j)-tol-a(j))/v(j),(hi(j)+tol-a(j))/v(j)]);
            tlo = max(tlo,ts(1)); thi = min(thi,ts(2));
            if tlo > thi, hit=false; return; end
        end
    end
    hit = true;
elseif D == 3 && numel(o) == 5
    r = o(4); h = o(5);
    if r < 0 || h < 0, error('rr:Obstacle','Cylinder radius/height must be nonnegative.'); end
    v = b-a; tol = 64*eps(max([1,abs(a),abs(b),abs(o)]));
    zlo = o(3)-h/2; zhi = o(3)+h/2; tlo = 0; thi = 1;
    if v(3) == 0
        if a(3) < zlo-tol || a(3) > zhi+tol, hit=false; return; end
    else
        ts=sort([(zlo-tol-a(3))/v(3),(zhi+tol-a(3))/v(3)]);
        tlo=max(0,ts(1)); thi=min(1,ts(2));
        if tlo>thi, hit=false; return; end
    end
    w = a(1:2)-o(1:2); vxy=v(1:2); denom=dot(vxy,vxy);
    if denom==0, t=tlo; else, t=min(thi,max(tlo,-dot(w,vxy)/denom)); end
    hit=norm(w+t*vxy)<=r+tol;
else
    error('rr:Obstacle','Unsupported obstacle representation for this dimension.');
end
end
