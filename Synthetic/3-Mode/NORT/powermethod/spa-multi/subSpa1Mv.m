function [ X ] = subSpa1Mv( subs, vals, V, tenSz )
% Spa * v1

nnZ = size(subs, 1);
X = zeros(tenSz(1), size(V,2));

for i = 1:nnZ
    sub1 = subs(i, 1);
    sub2 = subs(i, 2);
    sub3 = subs(i, 3);
    
    ai = vals(i);
    pi = sub2 + tenSz(2)*(sub3 - 1);
    
    ui = V(pi, :);
    
    X(sub1,:) = X(sub1,:) + ai*ui;
end

fprintf('should not be called\n');

end

