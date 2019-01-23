function [ vals ] = MakepartM2( subs, U, V, szTen )
% U2 V2

nnZ = size(subs, 1);

vals = zeros(nnZ, 1);

for i = 1:nnZ
    sub1 = subs(i, 1);
    sub2 = subs(i, 2);
    sub3 = subs(i, 3);
    
    u = U(sub2, :);
    v = V((sub1 - 1)*szTen(3) + sub3, :);
    
    vals(i) = u*v';
end

end

