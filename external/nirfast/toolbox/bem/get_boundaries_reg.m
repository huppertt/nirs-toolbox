function [inreg_a_1,inreg_b_1,inreg_c_1] = get_boundaries_reg(elem,allreg,reg);
%%finding all sub-boundaries in region 1
%% Subha 12/1/06
ee = find(elem(:,6) == (reg - 1));
ee1 = elem(ee,2:4);
ee1 = reshape(ee1',[],1);
inreg_a_1 = unique(ee1);
clear ee ee1;

inreg_b_1 = [];
inreg_c_1 = [];
for i = 1:length(allreg)
    ee = find(elem(:,6) == allreg(i));
    ee1 = elem(ee,2:4);
    ee1 = reshape(ee1',[],1);
    if (allreg(i) == 2)
        inreg_b_1 = unique(ee1);
    elseif (allreg(i) == 3)
        inreg_c_1 = unique(ee1);
    elseif (allreg(i) == 4)
        inreg_c_1 = unique(ee1);
    end
    clear ee ee1;
end
