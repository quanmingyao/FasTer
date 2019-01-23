function [ A ] = subSpa1u(Spa, U, sz)
% u1' * spa

n = sz(2);
A = zeros(size(U,2), sz(2)*sz(3));

for i = 1:length(Spa)
    B = U'*Spa{i};
    
    A(:, 1 + n*(i - 1):n*i) = B;
end

end