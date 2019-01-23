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
    lambda = [lambda, lambda, lambda];
end

Omega = (O ~= 0);
tenSz = size(O);
X0 = O;
X1 = X0;

c = 1;
obj = zeros(maxIter, 1);
RMSE = zeros(maxIter, 1);
Time = zeros(maxIter, 1);
Mused = zeros(maxIter, 1);

U = cell(3,1);
S = cell(3,1);
V = cell(3,1);
M = cell(3,1);

neg = 0;
for i = 1:maxIter
    tt = cputime;
    
    stepsize = 1;
    
    Y = X1 + (c - 1)/(c + 2)*(X1 - X0);
    
    G = Omega .* (Y - O);
    G = Y - (1/stepsize) * G;
    
    for m = 1:3
        [U{m}, S{m}, V{m}] = SVT_t(Unfold(G, tenSz, m), lambda(m)/stepsize);
        M{m} = U{m}*S{m}*V{m}';
    end
    
    Y = getX(M{1}, M{2}, M{3}, tenSz, lambda);
    
    obj(i) = getObj(O, Omega, Y, S{1}, S{2}, S{3}, lambda);
    
    if(i <= 1)
        delta = inf;
    else
        delta = (obj(i - 1) - obj(i))/obj(i);
    end
    
    if(delta < 0)
        c = 1;
        
        neg = neg + 1;
        X0 = X1;
    else
        c = c + 1;
        
        X0 = X1;
        X1 = Y;
    end
    
    fprintf('iter:%d, obj:(%.2d,%.2d,%2d), rnk:(%d,%d,%d) \n', ...
        i, obj(i), delta, c, nnz(S{1}), nnz(S{2}), nnz(S{3}));
    
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
    
    if(delta > 0 && delta < tol)
        break;
    end
    
    if(neg > 5)
        break;
    end
end

out.S = {diag(S{1}), diag(S{2}), diag(S{3})};
out.obj = obj(1:i);
out.RMSE = RMSE(1:i);
out.rank = [nnz(S{1}), nnz(S{2}), nnz(S{3})];
out.Time = cumsum(Time(1:i));
out.Mused = Mused(1:i);

end

%% ------------------------------------------------------------------------
function [ obj ] = getObj(O, Omega, X, S1, S2, S3, lambda)

obj = (O - X) .* Omega;
obj = (1/2) * sum(obj(:).^2);
obj = obj + lambda(1)*sum(S1(:));
obj = obj + lambda(2)*sum(S2(:));
obj = obj + lambda(3)*sum(S3(:));

end


