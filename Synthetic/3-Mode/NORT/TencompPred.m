function [ RMSE ] = TencompPred( U, S, V, wgt, test, tenSz )

for m = 1:length(U)
    U{m} = U{m} * diag(wgt(m)*S{m});
end

pred = subMakrPrat(test.subs, U, V, tenSz);
pred = pred - test.vals;

RMSE = sqrt(sum(pred.^2) / length(pred));

end

