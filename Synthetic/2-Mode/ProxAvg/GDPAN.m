function [ X1, out ] = GDPAN( O, lambda, theta, para, tstData )

if(isfield(para, 'maxIter'))
    maxIter = para.maxIter;
else
    maxIter = 5000;
end

if(isfield(para, 'tol'))
    tol = para.tol;
else
    tol = 1e-3;
end

if(isfield(para, 'regType'))
    regType = para.regType;
else
    regType = 2;
end

% if(isfield(para, 'maxR'))
%     maxR = para.maxR;
% else
%     maxR = [5, 5];
% end

if(length(lambda) == 1)
    lambda = [lambda, lambda];
end

if(length(theta) == 1)
    theta = [theta, theta];
end

if(isfield(para, 'X'))
    X1 = para.X;
else
    X1 = zeros(size(O));
end

Omega = (O ~= 0);
tenSz = size(O);

obj = zeros(maxIter, 1);
RMSE = zeros(maxIter, 1);
Time = zeros(maxIter, 1);
Mused = zeros(maxIter, 1);

M = cell(2, 1);
[U, S, V] = subInitUSV(tenSz, para.maxR);

for i = 1:maxIter
    tt = cputime;
    
    stepsize = 1;
    
    G = Omega .* (X1 - O);
    G = X1 - (1/stepsize) * G;
    
    for m = 1:2
        [U{m}, S{m}, V{m}] = powerMethod_gen(Unfold(G, tenSz, m), V{m}, 10);
        S{m} = proximalRegC_warpper(diag(S{m}), lambda(m)/stepsize, theta(m), regType);
        S{m} = diag(S{m});
        
%         [U{m}, S{m}, V{m}] = GSVT_t(Unfold(G, tenSz, m), lambda(m)/stepsize, ...
%             theta, regType);

%         S{m} = diag(S{m});
%         S{m}(min(para.maxR(m)+1,nnz(S{m})):end) = 0;
%         S{m} = diag(S{m});
        
        M{m} = U{m}*S{m}*V{m}';
    end
    
    X1 = getX(M, tenSz, lambda);
    
    obj(i) = getObj(O, Omega, X1, S, lambda, theta, regType);
    
    if(i <= 1)
        delta = inf;
    else
        delta = (obj(i - 1) - obj(i))/obj(i);
    end
    
    Time(i) =  cputime - tt;
    
    fprintf('iter:%d, obj:(%.2d,%.2d), rnk:(%d, %d) \n', ...
        i, obj(i), delta, nnz(S{1}), nnz(S{2}));
    temp = memory();
    Mused(i) = temp.MemUsedMATLAB/(1024^2);
    
    if(exist('tstData', 'var'))
        RMSEi = X1 - tstData;
        RMSEi = RMSEi.* (tstData ~= 0);
        RMSEi = sum(RMSEi(:).^2);
        RMSEi = sqrt(RMSEi/nnz(tstData));
        
        RMSE(i) = RMSEi;
        
        fprintf('testing:%.3d \n', RMSEi);
    end
    
    if(abs(delta) < tol)
        break;
    end
end

out.S = {diag(S{1}), diag(S{2})};
out.obj = obj(1:i);
out.RMSE = RMSE(1:i);
out.rank = [nnz(S{1}), nnz(S{2})];
out.Time = cumsum(Time(1:i));
out.Mused = Mused(1:i);

end

%% ------------------------------------------------------------------------
function [ obj ] = getObj(O, Omega, X, S, lambda, theta, regType)

obj = (O - X) .* Omega;
obj = (1/2) * sum(obj(:).^2);
obj = obj + funRegC(diag(S{1}), nnz(S{1}), lambda(1), theta, regType);
obj = obj + funRegC(diag(S{2}), nnz(S{2}), lambda(2), theta, regType);

end

%% ------------------------------------------------------------------------
function [ X ] = getX(M, tenSz, lambda)

alpha1 = lambda(1)/sum(lambda);
alpha2 = lambda(2)/sum(lambda);

X = alpha1 * Fold(M{1}, tenSz, 1) + alpha2 * Fold(M{2}, tenSz, 2);

end

%% ------------------------------------------------------------------------
function [U, S, V] = subInitUSV(tenSz, maxR)

U = cell(2, 1);
S = cell(2, 1);
V = cell(2, 1);

for i = 1:2
    U{i} = randn(tenSz(i), maxR(i));
    [U{i}, ~] = qr(U{i}, 0);
    
    S{i} = eye(maxR(i), maxR(i));
     
    V{i} = zeros(prod(tenSz)/tenSz(i), maxR(i));
end

end

