function [ x ] = subMulti3M1u( U3, V3, u, sz )
% u1 (U3 V3)

V = zeros(sz(2), size(V3,2));

for i = 1:size(U3, 2)
%     v2 = V3(:, i);
    v2 = reshape(V3(:, i), sz(1), sz(2));
    v2 = u'*v2;
    
    V(:,i) = v2';
end

x = U3*V';
x = reshape(x', 1, numel(x));

end

