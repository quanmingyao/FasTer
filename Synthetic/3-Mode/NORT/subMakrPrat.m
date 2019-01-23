function [ prti ] = subMakrPrat(subs, U0, V0, tenSz)

prti = 0;
prtm = zeros(size(subs,1), 1);

makeFunc = {@MakepartM1_c, @MakepartM2_c, @MakepartM3_c};
for m = 1:length(makeFunc)
    makeFunc{m}(subs, U0{m}, V0{m}, tenSz, prtm);
    prti = prti + prtm;
end

end