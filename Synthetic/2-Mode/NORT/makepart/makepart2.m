function [ part ] = makepart2( row, col, U2, S2, V2, sz )
% mode-2

n = sz(3);
part = cell(length(row), 1);

V = U2*diag(S2);

for i = 1:length(row)
    U = V2(i:n:size(V2,1), :);
    
    part{i} = sparse_inp(U', V', row{i}, col{i})';
end


end

