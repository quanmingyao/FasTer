function [ U0, S0, V0, out ] = Tencomp_acc( subs, vals, tenSz, lambda, ...
    theta, para )

if(isfield(para, 'maxIter'))
    maxIter = para.maxIter;
else
    maxIter = 5000;
end

if(isfield(para, 'maxR'))
    maxR = para.maxR;
else
    maxR = [5, 5, 5];
end

if(length(lambda) == 1)
    lambda = [lambda, lambda, lambda];
end

if(length(theta) == 1)
    theta = [theta, theta, theta];
end

[U0, S0, V0, wgt] = subInitUSV(tenSz, lambda);
U1 = U0;
V1 = V0;
S1 = S0;

cc = 1;
obj = zeros(maxIter, 1);
Time = zeros(maxIter, 1);
RMSE = zeros(maxIter, 1);
Mused = zeros(maxIter, 1);

% update sparse part
prt0 = subMakrPrat(subs, U0, V0, tenSz);

for i = 1:maxIter
    tt = cputime;
    
%     stepSz = sqrt(i)/(sqrt(i) + 3);
    
    stepSz = 1;
    
    beta0 = (cc - 1)/(cc + 2);
    beta1 = 1 + beta0;
    beta0 = - beta0;
    
    % merge S into U
    for m = 1:length(tenSz)
        U1{m} = U1{m}*diag(wgt(m)*S1{m});
    end
    
    % update sparse part
    prt1 = subMakrPrat(subs, U1, V1, tenSz);
    
    obj(i) = subGetObj(vals, prt1, S1, lambda, theta, para.regType);
    
    if(i == 1 || obj(i) < obj(i - 1))
        cc = cc + 1;
    else
        cc = 1;
    end
    
    prti = vals - (beta0*prt0 + beta1*prt1);
    prti = prti/stepSz;
    
    % SVT on each mode
    [Ut, Vt, R] = subScaleFactors(U1, V1, beta1, U0, V0, beta0, maxR);    
    [Ut, St, Vt] = subProxAvg(subs, prti, Ut, Vt, R, tenSz, maxR, ...
        lambda/stepSz, theta, para.regType);
    
    U0 = U1;
    S0 = S1;
    V0 = V1;
    U1 = Ut;
    S1 = St;
    V1 = Vt;
    prt0 = prt1;
    
    
    
    if(i == 1)
        delta = inf;
    else
        delta = abs(obj(i) - obj(i - 1))/obj(i);
    end
    
    fprintf('iter:%d,obj:(%.2d,%.1d),ac:%d,sz:%.1d,rnk:(%d,%d,%d)\n',...
        i, obj(i), delta, cc, stepSz, length(S1{1}), length(S1{2}), length(S1{3}));
    temp = memory();
    Mused(i) = temp.MemUsedMATLAB/(1024^2);
    
    Time(i) = cputime - tt;
    
    if(isfield(para, 'test'))
        RMSE(i) = TencompPred( U1, S1, V1, wgt, para.test, tenSz );
        fprintf('test: %.2d \n', RMSE(i));
    end
    
    
    if(i > 1 && delta < para.tol)
        break;
    end
end

out.obj = obj(1:i);
out.RMSE = RMSE(1:i);
out.Time = cumsum(Time(1:i));
out.S = S0;
out.Mused = Mused(1:i);

end

%% ------------------------------------------------------------------------
function [Ut, Vt, R] = subScaleFactors(U1, V1, beta1, U0, V0, beta0, maxR)

Ut = cell(length(U1), 1);
Vt = cell(length(V1), 1);
for m = 1:length(U1)
    Ut{m} = cat(2, beta1*U1{m}, beta0*U0{m});
    Vt{m} = cat(2, V1{m}, V0{m});
end

R = cell(length(U1), 1);

for m = 1:length(U1)
    [R{m}, ~] = qr(Ut{m}, 0);
    colM = min(maxR(m), size(R{m},2));
    R{m} = R{m}(:, 1:colM);
end

end

%% ------------------------------------------------------------------------
function [U1, S1, V1] = subProxAvg(subs, prt, Ut, Vt, R, ...
    tenSz, maxR, lambda, theta, regType)

pwPara.maxIter = 1;

U1 = cell(length(maxR), 1);
S1 = cell(length(maxR), 1);
V1 = cell(length(maxR), 1);

powerFunc = {@PowerMethod1, @PowerMethod2, @PowerMethod3};

for m = 1:length(tenSz)
    [ U1{m}, S1{m}, V1{m} ] = powerFunc{m}( subs, prt, Ut, Vt, R{m}, ...
        tenSz, pwPara );
    S1{m} = diag(S1{m});
    S1{m} = proximalRegC_warpper(S1{m}, lambda(m), theta(m), regType);

    [ U1{m}, S1{m}, V1{m}] = filterBase(U1{m}, S1{m}, V1{m}, maxR(m));
end

end

%% ------------------------------------------------------------------------
function [ obj ] = subGetObj(vals, pred, S, lambda, theta, regType)

obj = (1/2)*sum((vals - pred).^2);

for i = 1:length(S)
    obj = obj + funRegC(S{i}, length(S{i}), lambda(i), theta(i), regType);
end

end

%% ------------------------------------------------------------------------
function [U, S, V, alpha] = subInitUSV(tenSz, lambda)

U = cell(length(tenSz), 1);
S = cell(length(tenSz), 1);
V = cell(length(tenSz), 1);

for i = 1:length(tenSz)
    U{i} = randn(tenSz(i), 1);
    [U{i}, ~] = qr(U{i}, 0);
    
    S{i} = ones(1, 1);
     
    V{i} = zeros(prod(tenSz)/tenSz(i), 1);
end

alpha = lambda/sum(lambda);

end