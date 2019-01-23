function [ x ] = subMulti2M3v( U2, V2, v, sz )
% (U2 V2) v3

x = 0;
v = reshape(v, sz(1), sz(2));

for j = 1:size(V2, 2)
%     u1 = U2(:, j);
%     v1 = V2(:, j);
    v1 = reshape(V2(:, j), sz(3), sz(1));
    
    x = x + v1 * (v * U2(:, j));
end

end

