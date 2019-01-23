function [ tra, tst ] = SynData2M( tenSz, rnk, ratio)

I1= tenSz(1);
I2= tenSz(2);
I3= tenSz(3);

A = randn(I1,rnk);
% [A, ~] = qr(A, 0);
B = randn(I2,rnk);
% [B, ~] = qr(B, 0);
C = randn(I3,rnk);
% [C, ~] = qr(C, 0);
% s = prod(tenSz)*rand(rnk,1);

X = zeros(I1, I2, I3);

for i = 1:rnk
    x = kron(A(:,i)*B(:,i)', C(:,i));
    x = reshape(x, [I1,I2,I3]);
    
    X = X + x;
end

X = X + randn(size(X)) * 0.1;

tra = cell(I3, 1);
tst = struct;
for i = 1:I3
    temp = (rand(I1,I2) < ratio) ;
    
    tra{i} = temp .* X(:,:,i);
    tra{i} = sparse(tra{i});
    
    temp = (rand(I1,I2) < 0.005);
    temp = temp .* X(:,:,i);
    temp = sparse(temp);
    
    [tst.row{i}, tst.col{i}, tst.data{i}] = find(temp);
end

end

