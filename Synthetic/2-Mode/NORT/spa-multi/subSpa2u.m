function [ A ] = subSpa2u( Spa, U, sz )
% u2' * spa

A = zeros(size(U,2), sz(1)*sz(3));
n = size(A,2);

for i = 1:length(Spa)
    B = Spa{i}*U;
    
    A(:, i:length(Spa):n) = B';
end

end