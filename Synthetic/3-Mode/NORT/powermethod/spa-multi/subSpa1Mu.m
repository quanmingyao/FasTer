function [ X ] = subSpa1Mu( subs, vals, U, tenSz )
% u1' * Spa 

nnZ = size(subs, 1);
X = zeros(size(U,2), tenSz(2)*tenSz(3));

for i = 1:nnZ
    sub1 = subs(i, 1);
    sub2 = subs(i, 2);
    sub3 = subs(i, 3);
    
    ui = U(sub1, :);
    
    ai = vals(i);
    pi = sub2 + tenSz(2)*(sub3 - 1);
    
    X(:,pi) = X(:,pi) + ai*ui';
end

fprintf('should not be called\n');

end

