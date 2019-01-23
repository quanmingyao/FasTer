function [ RMSE ] = TencompPred( U, S, V, wgt, test, tenSz )

pred1 = makepart1( test.row, test.col, U{1}, S{1}*wgt(1), V{1}, tenSz );
pred2 = makepart2( test.row, test.col, U{2}, S{2}*wgt(2), V{2}, tenSz );

nnzR = 0;
RMSE = 0;
for i = 1:length(pred1)
    error = test.data{i} - pred1{i} - pred2{i};
    RMSE = RMSE + sum(error.^2);
    nnzR = nnzR + length(error);
end
RMSE = sqrt(RMSE/nnzR);

end

