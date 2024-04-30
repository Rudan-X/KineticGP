function [value,isterm,dir] = reaching_steadyA(t,~)
global Gs_VEL
if size(Gs_VEL,1)>2
    value      = abs((Gs_VEL(end,2)-Gs_VEL(end-1,2))/(Gs_VEL(end,1)-Gs_VEL(end-1,1)))<1e-4;
else
    value=false;
end

isterm = 1;
dir = 0;