function [ U1, S1, V1, out] = Tencomp_acc( data, lambda, theta, para )

if(isfield(para, 'maxR'))
    maxR = para.maxR;
else
    maxR = [5,5];
end

if(isfield(para, 'maxIter'))
    maxIter = para.maxIter;
else
    maxIter = 100;
end

if(length(lambda) == 1)
    lambda = [lambda, lambda];
end

wgt = [lambda(1)/sum(lambda), lambda(2)/sum(lambda)];
tenSz = [size(data{1},1), size(data{1},2), length(data)];

[U0, S0, V0, U1, S1, V1, row, col, obv, spa] = InitUSV(data);
Ui = U0;
Si = S0;
Vi = V0;

clear data;

pwPara.maxIter = 1;
obj = zeros(maxIter, 1);
Time = zeros(maxIter, 1);
RMSE = zeros(maxIter, 1);
Mused = zeros(maxIter, 1);

cc = 1;
prt0 = MakePart( row, col, U0, S0, V0, tenSz, wgt );
prt1 = MakePart( row, col, U1, S1, V1, tenSz, wgt );

for i = 1:maxIter
    tt = cputime;
    
    stepSz = sqrt(i)/(sqrt(i) + 1);
    % stepSz = 1;
    
    beta0 = (cc - 1)/(cc + 2);
    beta1 = 1 + beta0;
    beta0 = - beta0;
    
    % update sparse part
    prti = cell(length(prt0), 1);
    for p = 1:length(prt0)
        prti{p} = obv{p} - (beta1 * prt1{p} + beta0 * prt0{p});
    end
    spa = copyToSpa(prti, spa, stepSz);
    
    % low-rank part
    [Ut, Vt] = scaleFactor(U1, S1, V1, beta1, U0, S0, V0, beta0, wgt);
    powerFunc = {@PowerMethod1, @PowerMethod2};
    for m = 1:2
        R = warmStart(U1{m}, U0{m}, maxR);
        [Ui{m}, Si{m}, Vi{m}] = powerFunc{m}( spa, Ut{1}, Vt{1}, Ut{2}, ...
            Vt{2}, R, tenSz, pwPara);
        Si{m} = diag(Si{m});
        Si{m} = proximalRegC_warpper(Si{m}, lambda(m)/stepSz, theta, para.regType);
        [Ui{m}, Si{m}, Vi{m}] = filterBase(Ui{m}, Si{m}, Vi{m});
    end
    
    
    % make up sparse
    prti = MakePart( row, col, Ui, Si, Vi, tenSz, wgt );
    
    obj(i) = getObj(obv, prti, Si, lambda, theta);
    
    if(i == 1 || obj(i) < obj(i - 1))
        cc = cc + 1;
    else
        cc = 1;
    end
    
    U0 = U1;
    S0 = S1;
    V0 = V1;
    U1 = Ui;
    S1 = Si;
    V1 = Vi;

    prt0 = prt1;
    prt1 = prti;
    
    Time(i) = cputime - tt;
    temp = memory();
    Mused(i) = temp.MemUsedMATLAB/(1024^2);
    
    if(i == 1)
        delta = inf;
    else
        delta = abs(obj(i) - obj(i - 1))/obj(i);
    end
    
    fprintf('iter:%d, obj:(%.2d,%.2d),rnk:(%d,%d),ac:%d,sz::%.1d\n', ...
        i, obj(i), delta, length(Si{1}), length(Si{2}), cc, stepSz);
    
    if(isfield(para, 'test'))
        RMSE(i) = TencompPred( U1, S1, V1, wgt, para.test, tenSz );
        fprintf('RMSE:%.2d\n', RMSE(i));
    end
    
    if(delta < para.tol)
        break;
    end
end

S1{1} = wgt(1)*S1{1};
S1{2} = wgt(2)*S1{2};

out.obj = obj(1:i);
out.Time = cumsum(Time(1:i));
out.RMSE = RMSE(1:i);
out.rank = [nnz(S1{1}), nnz(S1{2})];
out.Mused = Mused(1:i);

end

%% ------------------------------------------------------------------------
function [Ut, Vt] = scaleFactor(U1, S1, V1, beta1, ...
    U0, S0, V0, beta0, wgt)

Ut = cell(length(U1), 1);
Vt = cell(length(V1), 1);

for m = 1:length(Ut)
    U1{m} = U1{m}*diag(beta1*wgt(m)*S1{m});
    U0{m} = U0{m}*diag(beta0*wgt(m)*S0{m});
    
    Ut{m} = cat(2, U1{m}, U0{m});
    Vt{m} = cat(2, V1{m}, V0{m});
end

end

%% ------------------------------------------------------------------------
function [U0, S0, V0, ...
    U1, S1, V1, ...
    row, col, obv, spa] = InitUSV( data )

sz = [size(data{1},1), size(data{1},2), length(data)];

U0{1} = randn(sz(1), 1);
[U0{1}, ~] = qr(U0{1}, 0);
S0{1} = zeros(1, 1);
V0{1} = zeros(sz(2)*sz(3),1);

U0{2} = randn(sz(2), 1);
[U0{2}, ~] = qr(U0{2}, 0);
S0{2} = zeros(1, 1);
V0{2} = zeros(sz(1)*sz(3), 1);

row = cell(length(data), 1);
col = cell(length(data), 1);
obv = cell(length(data), 1);
for i = 1:length(data)
    [row{i}, col{i}, obv{i}] = find(data{i});
end

spa = data;

U1 = U0;
V1 = V0;
S1 = S0;

end

%% ------------------------------------------------------------------------
function [ obj ] = getObj(obv, spa, S, lambda, theta)

obj = 0;
for i = 1:length(spa)
    obj = obj + 0.5*sum((obv{i} - spa{i}).^2);
end

obj = obj + funRegC(S{1}, length(S{1}), lambda(1), theta, 2);
obj = obj + funRegC(S{2}, length(S{2}), lambda(2), theta, 2);

end

%% ------------------------------------------------------------------------
function [U, S, V] = filterBase(U, S, V)

S = S(1:nnz(S));
U = U(:,1:length(S));
V = V(:,1:length(S));

end




