function [v, w] = orthonormal_basis(u)

u = u(:) / max(norm(u), 1e-12);

if abs(u(1)) < 0.9
    a = [1;0;0];
else
    a = [0;1;0];
end

v = cross(u, a); v = v / max(norm(v), 1e-12);
w = cross(u, v); w = w / max(norm(w), 1e-12);

v = v.';
w = w.';
end
