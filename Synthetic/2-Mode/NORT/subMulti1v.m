function [ A ] = subMulti1v( U1, V1, U2, V2, V, sz )
% (U1 V1' + U2 V2')*V1

A = zeros(size(U1,1), size(V,2));

for i = 1:size(V,2)
    v = V(:,i);
    
    % mode-2
    a = fastMulti2v( U2, V2, v, sz);
    % mode-1
    a = a + U1*(V1'*v);
    
    A(:,i) = a;
end

end

