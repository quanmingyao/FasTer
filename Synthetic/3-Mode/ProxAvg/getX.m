function [ X ] = getX(M1, M2, M3, tenSz, lambda)

alpha1 = lambda(1)/sum(lambda);
alpha2 = lambda(2)/sum(lambda);
alpha3 = lambda(3)/sum(lambda);

X = alpha1 * Fold(M1, tenSz, 1) ...
    + alpha2 * Fold(M2, tenSz, 2) ...
    + alpha3 * Fold(M3, tenSz, 3);

end