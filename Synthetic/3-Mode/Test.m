clear; clc;

tenSz = [10, 10, 10];
tenSz = 20*tenSz;

core_dims = 5*[1,1,1];
fraction = (tenSz(3)/core_dims(3)) * sum(tenSz)*log(prod(tenSz))/prod(tenSz);

for repeat = 1:1
    
maxNumCompThreads(1);
    
rng(repeat);   

[X, d] = synRand(tenSz, core_dims );

subs = makeOmegaSet(tenSz, round(fraction * prod(tenSz)));  
vals = getValsAtIndex_mex(subs', X.G, X.U{1}', ...
    X.U{2}', X.U{3}'); 
vals = vals + randn(size(vals))*0.01;

testsubs = makeOmegaSet(tenSz, round(0.1* fraction * prod(tenSz)));  
testvals = getValsAtIndex_mex(testsubs', X.G, X.U{1}', ...
    X.U{2}', X.U{3}'); 

X = ktensor(d, X.U);
X = full(X);
X = double(X);

para.maxIter = 5000;
para.tol = 1e-5;

gridLambda = 0.1 * 2.^(2:-1:-2);

for l = 1:length(gridLambda)
    para.maxR = core_dims;
    para.test.subs = testsubs;
    para.test.vals = testvals;
    
    memUsed = memory();
    memUsed = memUsed.MemUsedMATLAB/(1024^2);

    traX = zeros(size(X));
    traIdx = sub2ind(tenSz, subs(:,1), subs(:,2), subs(:,3));
    traX(traIdx) = vals;
    clear traIdx;

    tstX = zeros(size(X));
    tstIdx = sub2ind(tenSz, testsubs(:,1), testsubs(:,2), testsubs(:,3));
    tstX(tstIdx) = testvals;
    clear tstIdx;
    
    lambda = gridLambda(l);
  
    method = 1;
    para.tol = 1e-6;
    [~, out{method,l}] = ProxAvg(traX, lambda, para, tstX);
    
    para.maxR = [5,5,5];
    
    for regType = 1:1
        para.regType = regType;
        
        switch(regType)
            case 1
                lambda3 = 4*lambda;
                theta = 2*lambda3;
            case 2
                lambda3 = lambda;
                theta = sqrt(lambda3);
            case 3
                lambda3 = lambda/4;
                theta = 5;
        end
        
        method = method + 1;
        para.tol = 1e-5;
        [~, out{method,l}] = GDPAN(traX, lambda3, theta, para, tstX);

        method = method + 1;
        
        subs1 = subs;
        vals1 = vals;
        [~,~,~,out{method,l}] = sNORT( subs1, vals1, tenSz, lambda3, theta, para );
        
        clear subs1 vals1;

        method = method + 1;
        
        subs1 = subs;
        vals1 = vals;
        [~,~,~,out{method,l}] = NORT( subs1, vals1, tenSz, lambda3, theta, para );
    
        clear subs1 vals1;
    end
    
    clear lambda3;
       
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
    
    clear i j;
end


clear para subs traX tstX vals X;
clear lambda gridLambda method d testsubs testvals;

save( strcat(num2str(tenSz(1)), '-', num2str(repeat), '.mat') );

end



