function [ x ] = fastMulti1v( U1, V1, v, sz)
% (U1 V1) v2

x = 0;
v = reshape(v, sz(3), sz(1));

for j = 1:size(V1, 2)
%     u1 = U1(:, j);
%     v1 = V1(:, j);
    v1 = reshape(V1(:, j), sz(2), sz(3));
    
    x = x + v1*(v * U1(:, j));
end

end

