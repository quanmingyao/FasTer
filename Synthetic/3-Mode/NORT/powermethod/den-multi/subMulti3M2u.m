function [ x ] = subMulti3M2u( U3, V3, u, sz )
% u2 (U3 V3)

V = zeros(sz(1), size(V3,2));

for i = 1:size(U3, 2)
%     v2 = V3(:, i);
    v2 = reshape(V3(:, i), sz(1), sz(2));
    v2 = v2 * u;
    
    V(:,i) = v2;
end

x = U3*V';
x = reshape(x, 1, numel(x));

end

