function [ U, S, V, out] = Tencomp( data, lambda, theta, para )

if(isfield(para, 'maxR'))
    maxR = para.maxR;
else
    maxR = 5;
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

[U, S, V, row, col, obv, Spa] = InitUSV(data);
Ui = U;

clear data;

pwPara.maxIter = 3;
obj = zeros(maxIter, 1);
Time = zeros(maxIter, 1);
RMSE = zeros(maxIter, 1);
Mused = zeros(maxIter, 1);

for i = 1:maxIter
    tt = cputime;
    
    % makeup sparse part
    prt = MakePart( row, col, U, S, V, tenSz, wgt );
    for p = 1:length(prt)
        prt{p} = obv{p} - prt{p};
    end
    
    obj(i) = getObj(prt, S, lambda, theta);
    
    stepSz = 1;
    Spa = copyToSpa(prt, Spa, stepSz);
    
    U1 = U{1}*diag(wgt(1)*S{1});
    U2 = U{2}*diag(wgt(2)*S{2});
    V1 = V{1};
    V2 = V{2};
    
    % proximal-avg    
    R = warmStart(U1, Ui{1}, maxR);
    [ Ui{1}, S1, V{1} ] = PowerMethod1( Spa, U1, V1, U2, V2, R, ...
        tenSz, pwPara );
    S1 = diag(S1);
    S1 = proximalRegC_warpper(S1, lambda(1)/stepSz, theta, para.regType);
    [ U{1}, S{1}, V{1}, Ui{1} ] = filterBase(S1, V{1}, Ui{1});
     
    R = warmStart(U2, Ui{2}, maxR);
    [ Ui{2}, S2, V{2} ] = PowerMethod2( Spa, U1, V1, U2, V2, R, ...
        tenSz, pwPara );
    S2 = diag(S2);
    S2 = proximalRegC_warpper(S2, lambda(2)/stepSz, theta, para.regType);
    [ U{2}, S{2}, V{2}, Ui{2} ] = filterBase(S2, V{2}, Ui{2});
    
    Time(i) = cputime - tt;
    temp = memory();
    Mused(i) = temp.MemUsedMATLAB/(1024^2);
    if(i == 1)
        delta = inf;
    else
        delta = abs(obj(i) - obj(i - 1))/obj(i);
    end
    
    fprintf('iter:%d, obj:(%.2d,%.2d), rnk:(%d,%d)\n', ...
        i, obj(i), delta, length(S{1}), length(S{2}));
    if(isfield(para, 'test'))
        RMSE(i) = TencompPred( U, S, V, wgt, para.test, tenSz );
        fprintf('RMSE:%.2d\n', RMSE(i));
    end
    
    if(delta < para.tol)
        break;
    end
end

S{1} = wgt(1)*S{1};
S{2} = wgt(2)*S{2};

out.obj = obj(1:i);
out.Time = cumsum(Time(1:i));
out.RMSE = RMSE(1:i);
out.rank = [nnz(S{1}), nnz(S{2})];
out.Mused = Mused(1:i);


end

%% ------------------------------------------------------------------------
function [U, S, V, Ui] = filterBase(S, V, Ui)

Ui = Ui(:,1:nnz(S));
S = S(1:nnz(S));
U = Ui;
V = V(:,1:nnz(S));

end

%% ------------------------------------------------------------------------
function [ obj ] = getObj(prt, S, lambda, theta)

obj = 0;
for i = 1:length(prt)
    obj = obj + 0.5*sum((prt{i}).^2);
end

obj = obj + funRegC(S{1}, length(S{1}), lambda(1), theta, 2);
obj = obj + funRegC(S{2}, length(S{2}), lambda(2), theta, 2);

end


%% ------------------------------------------------------------------------
function [U, S, V, row, col, ...
    Obv, Spa] = InitUSV(data)

sz = [size(data{1},1), size(data{1},2), length(data)];

U{1} = randn(sz(1), 1);
[U{1}, ~] = qr(U{1}, 0);
S{1} = zeros(1, 1);
V{1} = zeros(sz(2)*sz(3), 1);

U{2} = randn(sz(2), 1);
[U{2}, ~] = qr(U{2}, 0);
S{2} = zeros(1, 1);
V{2} = zeros(sz(1)*sz(3), 1);

% prt = cell(length(data), 1);
row = cell(length(data), 1);
col = cell(length(data), 1);
Obv = cell(length(data), 1);
for i = 1:length(data)
    [row{i}, col{i}, Obv{i}] = find(data{i});
end

Spa = data;

end

%% ------------------------------------------------------------------------