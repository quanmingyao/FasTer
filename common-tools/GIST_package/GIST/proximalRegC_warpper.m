function [ s ] = proximalRegC_warpper( s, lambda, theta, regType )

switch(regType)
    case 1 % CAP
        s = proximalRegC(s, length(s), lambda, theta, 1);
    case 2 % Logrithm
        s = proximalRegC(s, length(s), lambda, theta, 2);
    case 3 % TNN
        z = max(s - lambda, 0);
        s(theta:end) = z(theta:end);
    otherwise
        assert(false);
end

end