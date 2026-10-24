function u = unit(v)

if ~isvector(v)
    error('输入必须为向量（行向量或列向量）');
end

vec_norm = norm(v);

if vec_norm < eps
    u = zeros(size(v));
    warning('输入向量为零向量，返回零向量');
else
    u = v / vec_norm;
end
end
