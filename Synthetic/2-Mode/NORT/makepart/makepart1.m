function [ part ] = makepart1( row, col, U1, S1, V1, sz )
% mode-1

n = sz(2);
part = cell(length(row), 1);

U = U1*diag(S1);

for i = 1:length(row)
    V = V1(1 + (i - 1)*n:i*n, :);
    
    part{i} = sparse_inp(U', V', row{i}, col{i})';
end


end

