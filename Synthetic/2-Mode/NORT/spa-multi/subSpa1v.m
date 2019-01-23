function [ A ] = subSpa1v(Spa, V, sz)
% spa * v1

n = size(Spa{1}, 2);
A = zeros(sz(1), size(V,2));

for i = 1:length(Spa)
    B = V(1 + (i - 1)*n : i*n, :);
    A = A + Spa{i}*B;
end

end