function [U, S, V] = powerMethod_gen(Z, R, maxIter)

Y = Z*R;
for i = 1:maxIter
    [Q, ~] = qr(Y, 0);
    Y = Z*(Z'*Q);
end

[U, S, V] = svd(Q'*Z, 'econ');
U = Q*U;

end