function [ U, S, V ] = SVT_t( Z, lambda, rnk )

Z = real(Z);

[U, S, V] = safeSVD( Z );

s = diag(S);
s = s - lambda;
s = s(s > 0);
svs = length(s);

U = U(:,1:svs);
V = V(:,1:svs);
S = diag(s);

end

%% ------------------------------------------------------------------------
function [U, S, V] = safeSVD( A )

try 
    [U, S, V] = svd(A, 'econ');
catch
    [U, S, V] = lansvd(A, floor(min(size(A))/4), 'L');
end

% [U, S, V] = svd(A, 'econ');

end