function [ U, S, V ] = filterBase( U, S, V, maxR )

nnzS = min(maxR, nnz(S));
S = S(1:nnzS);

U = U(:,1:nnzS);
V = V(:,1:nnzS);

end