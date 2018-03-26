function copytable2clip(tbl)
% This function copys a table object to the clipboard

if(ismember('model',tbl.Properties.VariableNames) && isa(tbl.model{1},'LinearModel'))
    tbl.model=[];
end


DD={};
flds = tbl.Properties.VariableNames;
for i=1:length(flds)
    DD{1,i}=flds{i};
    DD{2,i}=[char(13) '-----------------'];
    
    for j=1:height(tbl)
        if(iscell(tbl.(flds{i})))
            a=tbl.(flds{i}){j};
        else
            a=tbl.(flds{i})(j);
        end
        if(isnumeric(a))
            a=num2str(a);
        end
        if(islogical(a))
            if(a)
                a='TRUE';
            else
                a='FALSE';
            end
        end
        DD{2+j,i}=a;
    end
end


S=[];
for idx=1:size(DD,1)
    for idx2=1:size(DD,2)
        S=[S DD{idx,idx2}];
        if(idx2~=size(DD,2))
            S=[S char(9)];
        else
            S=[S char(10)];
        end
    end
    
end

clipboard('copy',S);