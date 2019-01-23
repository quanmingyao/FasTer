function [ part ] = MakePart( row, col, U, S, V, sz, wgt )

part1 = makepart1( row, col, U{1}, wgt(1)*S{1}, V{1}, sz );
part2 = makepart2( row, col, U{2}, wgt(2)*S{2}, V{2}, sz );

part = cell(length(row), 1);
for i = 1:length(row)
    part{i} = part1{i} + part2{i};
end

end

