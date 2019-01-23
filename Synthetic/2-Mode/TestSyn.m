clear; clc;
close all;

for repeat = 1:1

maxNumCompThreads(1);

rng(repeat);

tenSz = [5, 5];
tenSz = [50*tenSz, 5];

core_dims = [5,5,tenSz(3)];
fraction = core_dims(3) * sum(tenSz) * log(prod(tenSz)) / prod(tenSz);

A = randn(tenSz(1),tenSz(3));
% [A, ~] = qr(A, 0);
A = bsxfun(@minus, A, mean(A,1));
B = randn(tenSz(2),tenSz(3));
% [B, ~] = qr(B, 0);
B = bsxfun(@minus, B, mean(B,1));
C = randn(tenSz(3),tenSz(3));
% [C, ~] = qr(C, 0);
C = bsxfun(@minus, C, mean(C,1));


X = ktensor(rand(tenSz(3),1), A, B, C);
X = full(X);
X = double(X);
X = X + randn(size(X)) * 0.01;

clear tst A B C;

gridLambda = 2.^(0:-1:-4);

for l = 1:length(gridLambda)
    addpath(genpath('Tencomp'));
    
    traData = zeros(size(X));
    tstData = zeros(size(X));
    
    for i = 1:tenSz(3)
        temp = (rand(tenSz(1),tenSz(2)) < fraction) ;
        temp = temp .* X(:,:,i);

        traData(:,:,i) = temp;

        temp = (rand(tenSz(1),tenSz(2)) < 0.01);
        temp = temp .* X(:,:,i);
        
        tstData(:,:,i) = temp;
        
        temp = sparse(temp);

        [para.test.row{i}, para.test.col{i}, para.test.data{i}] = find(temp);
    end

    para.maxIter = 5000;

    lambda = gridLambda(l);
    
    %% --------------------------------------------------------------------
    method = 1;
    para.tol = 1e-5;
    temp = traData + 0.0;
    [ ~, out{method,l} ] = ProxAvg( temp, lambda, para, tstData );
    
    clear temp;
    
    %% --------------------------------------------------------------------
    for regType = 1:1
        para.regType = regType;
        
        switch(regType)
            case 1
                lambda2 = 4*lambda;
                theta = 2*lambda2;
            case 2
                lambda2 = 8*lambda;
                theta = sqrt(lambda2);
            case 3
                lambda2 = lambda/4;
                theta = 5;
        end
        
        method = method + 1;
        para.tol = 1e-5;
        para.maxR = [5, 5];

        [ ~, out{method,l} ] = GDPAN( traData, lambda2, theta, para, tstData );

        %% --------------------------------------------------------------------
        method = method + 1;
        
        traSp = cell(tenSz(3), 1);
        for i = 1:tenSz(3)
            traSp{i} = sparse(traData(:,:,i));
        end
        
        para.tol = 1e-5;

        [ ~, ~, ~, out{method,l}] = sNORT( traSp, lambda2, theta, para );

        clear traSp;
        
        %% --------------------------------------------------------------------
        method = method + 1;
        
        traSp = cell(tenSz(3), 1);
        for i = 1:tenSz(3)
            traSp{i} = sparse(traData(:,:,i));
        end
        
        para.tol = 1e-5;
        [ ~, ~, ~, out{method,l}] = NORT( traSp, lambda2, theta, para );
    end
    
    clear lambda2;
    
    RMSE = zeros(size(out));
    Time = zeros(size(out));
    for i = 1:size(out,1)
        for j = 1:size(out, 2)
            if(~isempty(out{i,j}))
                RMSE(i,j) = out{i,j}.RMSE(end);
                Time(i,j) = out{i,j}.Time(end);
            end
        end
    end
    
    clear i j ans regType theta;
    
    save(strcat(num2str(tenSz(1)),'-',num2str(repeat),'.mat'), 'RMSE', 'Time', 'out');
end

clear i l lambda temp para traData traSp tstData X gridLambda;

save( strcat(num2str(tenSz(1)),'-',num2str(repeat),'.mat') );

end




