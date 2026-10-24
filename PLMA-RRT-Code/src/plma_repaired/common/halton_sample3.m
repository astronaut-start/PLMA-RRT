function q = halton_sample3(k, space)

b = [2, 3, 5];
h = [radical_inverse(k, b(1)), ...
     radical_inverse(k, b(2)), ...
     radical_inverse(k, b(3))];

q = zeros(1,3);
q(1) = space(1,1) + h(1) * (space(1,2) - space(1,1));
q(2) = space(2,1) + h(2) * (space(2,2) - space(2,1));
q(3) = space(3,1) + h(3) * (space(3,2) - space(3,1));
end

function x = radical_inverse(i, base)

x = 0;
f = 1 / base;

while i > 0
    digit = mod(i, base);
    x = x + digit * f;
    i = floor(i / base);
    f = f / base;
end
end
