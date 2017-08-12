function [inreg_1,inreg_2,inreg_3,inreg_4] = get_region_nodes(melements,nregions)
%%calculates nodes in each region (eg: fat/glandular) using material
%%numbers in element list
%%for 3D element list
%%subha 8/17/06
%%currently programmed for maximum 4 regions
inreg_2 = [];
inreg_3 = [];
inreg_4 = [];

for i = 1:nregions
    ee = find(melements(:,5) == i);
    
    ee1 = melements(ee,2:4);
    ee1 = reshape(ee1',[],1);
    temp = unique(ee1);
    clear ee ee1;
    
    if (i==1)
        inreg_1 = temp;
    elseif (i == 2)
        inreg_2 = temp;
    elseif (i == 3)
        inreg_3 = temp;
    elseif (i == 4)
        inreg_4 = temp;
    end
    clear temp
end

