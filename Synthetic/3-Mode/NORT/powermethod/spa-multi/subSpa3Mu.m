function [ X ] = subSpa3Mu( subs, vals, U, tenSz )
% u3' * Spa

nnZ = size(subs, 1);
X = zeros(size(U,2), tenSz(1)*tenSz(2));

for i = 1:nnZ
    sub1 = subs(i, 1);
    sub2 = subs(i, 2);
    sub3 = subs(i, 3);
    
    ui = U(sub3, :);
    
    ai = vals(i);
    pi = sub1 + tenSz(1)*(sub2 - 1);

    X(:,pi) = X(:,pi) + ai*ui';
end

fprintf('should not be called\n');

end

