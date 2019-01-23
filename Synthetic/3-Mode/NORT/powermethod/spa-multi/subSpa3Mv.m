function [ X ] = subSpa3Mv( subs, vals, V, tenSz )
% Spa * v3

nnZ = size(subs, 1);
X = zeros(tenSz(3), size(V,2));

for i = 1:nnZ
    sub1 = subs(i, 1);
    sub2 = subs(i, 2);
    sub3 = subs(i, 3);
    
    ai = vals(i);
    pi = sub1 + tenSz(1)*(sub2 - 1);
    
    ui = V(pi, :);
    
    X(sub3,:) = X(sub3,:) + ai*ui;
end

fprintf('should not be called\n');

end

