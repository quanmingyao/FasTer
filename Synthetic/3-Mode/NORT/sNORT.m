function [ U0, S0, V0, out ] = Tencomp( subs, vals, tenSz, lambda, theta, para )

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

regType = para.regType;
[U0, S0, V0, wgt] = subInitUSV(tenSz, maxR, lambda);

Ui = U0;
prti = zeros(size(vals));

obj = zeros(maxIter, 1);
Time = zeros(maxIter, 1);
RMSE = zeros(maxIter, 1);
Mused = zeros(maxIter, 1);

for i = 1:maxIter
    tt = cputime;
    
    stepSz = 1;
    
    % merge S into U
    for m = 1:length(tenSz)
        U0{m} = U0{m}*diag(wgt(m)*S0{m});
    end
    
    % update sparse part
    prt = vals;
    makeFunc = {@MakepartM1_c, @MakepartM2_c, @MakepartM3_c};
    for m = 1:length(makeFunc)
        makeFunc{m}(subs, U0{m}, V0{m}, tenSz, prti);
        prt = prt - prti;
    end
    
    obj(i) = subGetObj(prt, S0, lambda, theta, regType);
    
    prt = prt/stepSz;
    
    % SVT on each mode
    [U0, S0, V0, Ui] = subProxAvg(subs, prt, U0, V0, Ui, tenSz, maxR, ...
        lambda/stepSz, theta, regType);
    
    if(i == 1)
        delta = inf;
    else
        delta = abs(obj(i) - obj(i - 1))/obj(i);
    end
    
    fprintf('iter:%d,obj:(%.2d,%.2d),rnk:(%d,%d,%d)\n',...
        i, obj(i), delta, length(S0{1}), length(S0{2}), length(S0{3}));
    
    Time(i) = cputime - tt;
    temp = memory();
    Mused(i) = temp.MemUsedMATLAB/(1024^2);
    
    if(isfield(para, 'test'))
        RMSE(i) = TencompPred( U0, S0, V0, wgt, para.test, tenSz );
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
function [U1, S1, V1, U0] = subProxAvg(subs, prt, U0, V0, Ui, tenSz, maxR, ...
    lambda, theta, regType)

pwPara.maxIter = 3;

U1 = cell(length(maxR), 1);
S1 = cell(length(maxR), 1);
V1 = cell(length(maxR), 1);

powerFunc = {@PowerMethod1, @PowerMethod2, @PowerMethod3};

for m = 1:length(tenSz)
    R = subWarmStart(U0{m}, Ui{m}, maxR(m));
    [ U1{m}, S1{m}, V1{m} ] = powerFunc{m}( subs, prt, U0, V0, R, ...
        tenSz, pwPara );
    S1{m} = diag(S1{m});
    S1{m} = proximalRegC_warpper(S1{m}, lambda(m), theta, regType);

    [ U1{m}, S1{m}, V1{m}] = filterBase(U1{m}, S1{m}, V1{m}, maxR(m));
end

end

%% ------------------------------------------------------------------------
function [ R ] = subWarmStart(U1, U0, maxR)

[R, ~] = qr([U1, U0], 0);
R = R(:, 1:min(maxR, size(R,2)));

end

%% ------------------------------------------------------------------------
function [ obj ] = subGetObj(val, S, lambda, theta, regType)

obj = (1/2)*sum(val.^2);

for i = 1:length(S)
    obj = obj + funRegC(S{i}, length(S{i}), lambda(i), theta, regType);
end

end

%% ------------------------------------------------------------------------
function [U, S, V, alpha] = subInitUSV(tenSz, maxR, lambda)

U = cell(length(tenSz), 1);
S = cell(length(tenSz), 1);
V = cell(length(tenSz), 1);

for i = 1:length(tenSz)
    U{i} = randn(tenSz(i), maxR(i));
    [U{i}, ~] = qr(U{i}, 0);
    
    S{i} = ones(maxR(i), 1);
     
    V{i} = zeros(prod(tenSz)/tenSz(i), maxR(i));
end

alpha = lambda/sum(lambda);

end