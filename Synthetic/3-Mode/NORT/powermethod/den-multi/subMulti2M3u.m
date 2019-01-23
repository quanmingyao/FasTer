function [ x ] = subMulti2M3u( U2, V2, u, sz )
% u3 (U2 V2)

V = zeros(sz(1), size(V2,2));

for i = 1:size(U2, 2)
%     v2 = V2(:, i);
    v2 = reshape(V2(:, i), sz(3), sz(1));
    v2 = u'*v2;
    
    V(:,i) = v2';
end

x = U2*V';
x = reshape(x', 1, numel(x));

end

