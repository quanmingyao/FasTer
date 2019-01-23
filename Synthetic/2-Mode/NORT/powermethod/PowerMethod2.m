function [ U, S, V ] = PowerMethod2( Spa, U1, V1, U2, V2, R, sz, para )
% 2nd mode

maxIter = para.maxIter;

Y = subMultiu(Spa, U1, V1, U2, V2, R, sz)';

i = 1;
while(1)
    [Q, ~] = qr(Y, 0);
    
    if(i == maxIter)
        break;
    else
        i = i + 1;
    end
    
    Y = subMultiv(Spa, U1, V1, U2, V2, Q, sz);
    Y = subMultiu(Spa, U1, V1, U2, V2, Y, sz)';
end

[U, S, V] = reduceSVD(Spa, U1, V1, U2, V2, Q, sz);

end

%% ------------------------------------------------------------------------
function [ A ] = subMultiu(Spa, U1, V1, U2, V2, U, sz)

% low-rank
A = subMulti2u( U1, V1, U2, V2, U, sz );

% sparse
A = A + subSpa2u( Spa, U, sz );

end

%% ------------------------------------------------------------------------
function [ A ] = subMultiv(Spa, U1, V1, U2, V2, V, sz)

% low-rank
A = subMulti2v( U1, V1, U2, V2, V, sz );

% sparse
A = A + subSpa2v(Spa, V, sz);

end

%% ------------------------------------------------------------------------
function [U, S, V] = reduceSVD(Spa, U1, V1, U2, V2, V, sz)

U = subMultiv(Spa, U1, V1, U2, V2, V, sz);
[U, R] = qr(U, 0);
[Ur, S, Vr] = svd(R, 'econ');
U = U*Ur;
V = V*Vr;

end