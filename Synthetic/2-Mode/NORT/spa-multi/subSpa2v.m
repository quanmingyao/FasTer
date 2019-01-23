function [ A ] = subSpa2v(Spa, V, sz)
% spa * v2

n = sz(3);
A = zeros(sz(2), size(V,2));

for i = 1:length(Spa)
    B = V(i:n:size(V,1), :);
    A = A + Spa{i}'*B;
end

end