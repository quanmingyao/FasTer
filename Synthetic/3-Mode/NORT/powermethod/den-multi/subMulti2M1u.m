function [ x ] = subMulti2M1u( U2, V2, u, sz )
% u1 (U2 V2)

V = zeros(size(V2,2), sz(3));

for i = 1:size(U2, 2)
%     v2 = V2(:, i);
    v2 = reshape(V2(:, i), sz(3), sz(1));
    v2 = v2*u;
    
    V(i, :) = v2';
end

x = U2*V;
x = reshape(x, 1, numel(x));

end

