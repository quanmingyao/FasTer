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

if(1 == length(lambda))
    lambda = [lambda, lambda, lambda];
end

if(length(theta) == 1)
    theta = [theta, theta, theta];
end

if(isfield(para, 'X'))
    X0 = para.X;
else
    X0 = zeros(size(O));
end

Omega = (O ~= 0);
tenSz = size(O);
X1 = X0;

obj = zeros(maxIter, 1);
RMSE = zeros(maxIter, 1);
Time = zeros(maxIter, 1);
Mused = zeros(maxIter, 1);

[U, S, V] = subInitUSV(tenSz, para.maxR);

for i = 1:maxIter
    tt = cputime;
    
    stepsize = 1;
    
    G = Omega .* (X1 - O);
    G = X1 - (1/stepsize) * G;
    
    S = cell(3, 1);
    M = cell(3, 1);
    
    for m = 1:3
        [U{m}, S{m}, V{m}] = powerMethod_gen(Unfold(G, tenSz, m), V{m}, 10);
        S{m} = proximalRegC_warpper(diag(S{m}), lambda(m)/stepsize, theta(m), para.regType);
        S{m} = diag(S{m});
        
%         [U{m}, S{m}, V{m}] = GSVT_t(Unfold(G, tenSz, m), ...
%             lambda(m)/stepsize, theta, para.regType);
%         S{m} = diag(S{m});
        
        M{m} = U{m}*S{m}*V{m}';
    end    
    
    X1 = getX(M{1}, M{2}, M{3}, tenSz, lambda);
    
    obj(i) = proxAvgObj(O, Omega, X1, S, lambda, theta, para.regType);
    
    if(i <= 1)
        delta = inf;
    else
        delta = (obj(i - 1) - obj(i))/obj(i);
    end
    
    fprintf('iter:%d, obj:(%.2d,%.2d) \n', i, obj(i), delta);
    
    Time(i) = cputime - tt;
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

out.S = {diag(S{1}), diag(S{2}), diag(S{3})};
out.obj = obj(1:i);
out.RMSE = RMSE(1:i);
out.rank = [nnz(S{1}), nnz(S{2}), nnz(S{3})];
out.Time = cumsum(Time(1:i));
out.Mused = Mused(1:i);

end


%% ------------------------------------------------------------------------
function [U, S, V] = subInitUSV(tenSz, maxR)

U = cell(length(tenSz), 1);
S = cell(length(tenSz), 1);
V = cell(length(tenSz), 1);

for i = 1:length(tenSz)
    U{i} = randn(tenSz(i), maxR(i));
    [U{i}, ~] = qr(U{i}, 0);
    
    S{i} = ones(maxR(i), maxR(i));
     
    V{i} = zeros(prod(tenSz)/tenSz(i), maxR(i));
end

end
