function [ x ] = subMulti3M1v( U3, V3, v, sz )
% (U3 V3) v1

x = 0;
v = reshape(v, sz(2), sz(3));

for j = 1:size(V3, 2)
%     u1 = U3(:, j);
%     v1 = V3(:, j);
    v1 = reshape(V3(:, j), sz(1), sz(2));
    
    x = x + v1 * (v * U3(:, j));
end

end

