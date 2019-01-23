function [ obj ] = proxAvgObj(O, Omega, X, S, lambda, theta, regType)

obj = (O - X) .* Omega;
obj = (1/2) * sum(obj(:).^2);

for m = 1:length(S)
    obj = obj + funRegC(S{m}, length(S{m}), lambda(m), theta, regType);
end

end