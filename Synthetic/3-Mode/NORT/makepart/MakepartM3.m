function [ vals ] = MakepartM3( subs, U, V, szTen )
% U3 V3

nnZ = size(subs, 1);

vals = zeros(nnZ, 1);

for i = 1:nnZ
    sub1 = subs(i, 1);
    sub2 = subs(i, 2);
    sub3 = subs(i, 3);
    
    u = U(sub3, :);
    v = V((sub2 - 1)*szTen(1) + sub1, :);
    
    vals(i) = u*v';
end

end

