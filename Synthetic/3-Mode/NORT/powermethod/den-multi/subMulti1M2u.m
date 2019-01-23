function [ x ] = subMulti1M2u( U1, V1, u, sz )
% u2 (U1 V1)

V = zeros(size(V1,2), sz(3));

for i = 1:size(V1, 2)
%     v1 = V1(:,i);
    v1 = reshape(V1(:,i), sz(2), sz(3));
    v1 = u'*v1;
    
    V(i, :) = v1;
end

x = U1*V;
x = reshape(x', 1, numel(x));

end

