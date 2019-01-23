function [ A ] = subMulti2u( U1, V1, U2, V2, U, sz )
% U2'*(U1 V1' + U2 V2')

A = zeros(size(U,2), size(V2,1));

for i = 1:size(U,2)
    u = U(:,i);
    
    % mode-1
    a = fastMulti1u( U1, V1, u, sz);
    % mode-2
    a = a + (u'*U2)*V2';
    
    A(i,:) = a;
end

end

