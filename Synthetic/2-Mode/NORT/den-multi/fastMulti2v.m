function [ x ] = fastMulti2v( U2, V2, v, sz)
% u1 (U2 V2)

x = 0;
v = reshape(v, sz(2), sz(3));

for j = 1:size(V2, 2)
%     u1 = U2(:, j);
%     v1 = V2(:, j);
    v1 = reshape(V2(:, j), sz(3), sz(1));
    
    x = x + v1'*(v' * U2(:, j));
end

end

