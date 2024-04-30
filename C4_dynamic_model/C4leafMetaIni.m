function [LeafMs,Vel,KVal] = C4leafMetaIni(PARAMS,KVALUES,VMAX)

intercept=PARAMS(2);
leafInis = zeros(1,7);
leafInis= LeafIni(intercept);
MetaInis = zeros(1,87);
[MetaInis,Vel,KVal]= C4Ini(PARAMS,KVALUES,VMAX);

for m = 1:7
    LeafMs(m) = leafInis(m);
end

for m=1:87
    LeafMs(7+m)= MetaInis(m);
end