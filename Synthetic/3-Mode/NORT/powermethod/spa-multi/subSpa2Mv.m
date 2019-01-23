function [ X ] = subSpa2Mv( subs, vals, V, tenSz )
% Spa * v2

nnZ = size(subs, 1);
X = zeros(tenSz(2), size(V,2));

for i = 1:nnZ
    sub1 = subs(i, 1);
    sub2 = subs(i, 2);
    sub3 = subs(i, 3);
    
    ai = vals(i);
    pi = sub3 + tenSz(3)*(sub1 - 1);
    
    ui = V(pi, :);
    
    X(sub2,:) = X(sub2,:) + ai*ui;
end

fprintf('should not be called\n');

end

