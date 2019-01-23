function [Spa] = copyToSpa(part, Spa, stepSz)

for i = 1:length(part)
    sparse_update(Spa{i}, part{i}/stepSz);
end

end