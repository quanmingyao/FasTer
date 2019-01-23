function [ A ] = subMulti1u( U1, V1, U2, V2, U, sz )
% U1'*(U1 V1' + U2 V2')

A = zeros(size(U,2), size(V1,1));

for i = 1:size(U,2)
    u = U(:,i);
    
    % mode-2
    a = fastMulti2u( U2, V2, u, sz);
    % mode-1
    a = a + (u'*U1)*V1';
    
    A(i,:) = a;
end

end

