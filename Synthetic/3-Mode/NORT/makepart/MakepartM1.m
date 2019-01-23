function [ vals ] = MakepartM1( subs, U, V, szTen )
% U1 V1

nnZ = size(subs, 1);

vals = zeros(nnZ, 1);

for i = 1:nnZ
    sub1 = subs(i, 1);
    sub2 = subs(i, 2);
    sub3 = subs(i, 3);
    
    u = U(sub1, :);
    v = V((sub3 - 1)*szTen(2) + sub2, :);
    
    vals(i) = u*v';
end

end

