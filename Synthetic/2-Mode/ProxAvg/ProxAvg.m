function [ X1, out ] = ProxAvg( O, lambda, para, tstData )

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

if(1 == length(lambda))
    lambda = [lambda, lambda];
end

Omega = (O ~= 0);
tenSz = size(O);
X0 = zeros(tenSz);
X1 = X0;

c = 1;
obj = zeros(maxIter, 1);
RMSE = zeros(maxIter, 1);
Time = zeros(maxIter, 1);
Mused = zeros(maxIter, 1);

stepsize = 1;
neg = 0;

for i = 1:maxIter
    tt = cputime;
    
    Y = X1 + (c - 1)/(c + 2)*(X1 - X0);
    
    G = Omega .* (Y - O);
    G = Y - (1/stepsize) * G;
    
    [U1, S1, V1] = SVT_t(Unfold(G, tenSz, 1), lambda(1)/stepsize);
    M1 = U1*S1*V1';
    
    [U2, S2, V2] = SVT_t(Unfold(G, tenSz, 2), lambda(2)/stepsize);
    M2 = U2*S2*V2';
    
    Y = getX(M1, M2, tenSz, lambda);
    
    obj(i) = getObj(O, Omega, Y, S1, S2, lambda);
    
    if(i <= 1)
        delta = inf;
    else
        delta = (obj(i - 1) - obj(i))/obj(i);
    end
    
    if(delta < 0)
        neg = neg + 1;
        c = 1;
        
        X0 = X1;
    else
        c = c + 1;
        
        X0 = X1;
        X1 = Y;
    end
    
    fprintf('iter:%d, obj:(%.2d,%.2d,%2d), rnk:(%d, %d) \n', ...
        i, obj(i), delta, c, nnz(S1), nnz(S2));
    
    Time(i) = cputime - tt;
    temp = memory();
    Mused(i) = temp.MemUsedMATLAB/(1024^2);
    if(exist('tstData', 'var'))
        RMSEi = X0 - tstData;
        RMSEi = RMSEi.* (tstData ~= 0);
        RMSEi = sum(RMSEi(:).^2);
        RMSEi = sqrt(RMSEi/nnz(tstData));
        
        RMSE(i) = RMSEi;
        
        fprintf('testing:%.3d \n', RMSEi);
    end
    
    if(delta > 0 && delta <= tol)
        break;
    end
    
    if(neg > 5)
        break;
    end
end

out.S = {S1, S2};
out.obj = obj(1:i);
out.RMSE = RMSE(1:i);
out.rank = [nnz(S1), nnz(S2)];
out.Time = cumsum(Time(1:i));
out.Mused = Mused(1:i);

end

%% ------------------------------------------------------------------------
function [ obj ] = getObj(O, Omega, X, S1, S2, lambda)

obj = (O - X) .* Omega;
obj = (1/2) * sum(obj(:).^2);
obj = obj + lambda(1)*sum(S1(:));
obj = obj + lambda(2)*sum(S2(:));

end

%% ------------------------------------------------------------------------
function [ X ] = getX(M1, M2, tenSz, lambda)

alpha1 = lambda(1)/sum(lambda);
alpha2 = lambda(2)/sum(lambda);

X = alpha1 * Fold(M1, tenSz, 1) + alpha2 * Fold(M2, tenSz, 2);

end

