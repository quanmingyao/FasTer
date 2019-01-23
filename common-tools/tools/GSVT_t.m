function [ U, S, V ] = GSVT_t( Z, lambda, theta, regType, rnk )
%% ------------------------------------------------------------------------
% exact solve low rank proximal step
% (1/2)*|X - Z|_F^2 + lambda |X|_theta
%% ------------------------------------------------------------------------
%  regtype = 1: Capped L1 regularizer 
%            2: Log Sum Penalty
%            3: TNN
%% ------------------------------------------------------------------------

if(exist('rnk', 'var'))
    rnk = min([rnk, min(size(Z))]);
    [U, S, V] = lansvd(Z, rnk, 'L');
else
    [m, n] = size(Z);
    if(m > n)
        [U, S, V] = svd(Z, 'econ');
    else
        [V, S, U] = svd(Z', 'econ');
    end
end

s = diag(S);
s = proximalRegC_warpper(s, lambda, theta, regType);
S = diag(s);

% switch(regType)
%     case 1 % CAP
%         s = proximalRegC(s, length(s), lambda, theta, 1);
%     case 2 % Logrithm
%         s = proximalRegC(s, length(s), lambda, theta, 2);
%     case 3 % TNN
%         z = max(s - lambda, 0);
%         s(theta:end) = z(theta:end);
%     otherwise
%         assert(false);
% end


% svs = sum(s > 1e-10);

% U = U(:,1:svs);
% V = V(:,1:svs);
% S = diag(s(1:svs));



end