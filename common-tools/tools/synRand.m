function [ X, d ] = synRand( n, k )

U1 = rand( n(1), k(1) );
U1 = bsxfun(@minus, U1, mean(U1, 1));
U2 = rand( n(2), k(2) );
U2 = bsxfun(@minus, U2, mean(U2, 1));
U3 = rand( n(3), k(3) );
U3 = bsxfun(@minus, U3, mean(U3, 1));

d = randn(k(1), 1);
d = d - mean(d);

C = zeros(k); 
for i = 1:length(d)
    C(i,i,i) = d(i);
end

X.U{1} = U1;
X.U{2} = U2;
X.U{3} = U3;
X.G = C;
    
end
