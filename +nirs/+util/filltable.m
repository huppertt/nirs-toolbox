function newtbl=filltable(newtbl,template,fillval)
% This function fills in the "new table" fields using the template table

lst=find(~ismember(template.Properties.VariableNames,newtbl.Properties.VariableNames));

for idx=1:length(lst)
    val=template.(template.Properties.VariableNames{lst(idx)})(1);
    if(iscell(val));
        val={fillval};
    else
        val=fillval;
    end;
    
    val=repmat(val,height(newtbl),1);
    newtbl=[newtbl table(val,'VariableNames',{template.Properties.VariableNames{lst(idx)}})];

end
return
