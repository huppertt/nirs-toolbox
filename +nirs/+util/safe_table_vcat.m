function out=safe_table_vcat(tbls,tbl2)

if(isempty(tbls))
    out=tbl2;
    return
end

if(~iscell(tbls) & nargin==2)
    tbls={tbls; tbl2};
end


flds={};
types=[];
for i=1:length(tbls)
    flds={flds{:} tbls{i}.Properties.VariableNames{:}};
    if(~isfield(tbls{i}.Properties,'VariableTypes'))
        typ={};
        for ii=1:length(tbls{i}.Properties.VariableNames)
            typ{ii}=class(tbls{i}.(tbls{i}.Properties.VariableNames{ii}));
        end
        types=[types(:); typ(:)];
    else
        types=[types(:); tbls{i}.Properties.VariableTypes(:)];
    end
end
[flds,~,id]=unique(flds);

types_safe={};
for i=1:length(flds)
    if(all(ismember(types(find(id==i)),types(min(find(id==i))))))
        types_safe{i}=types{min(find(id==i))};
    else
        types_safe{i}="mixed";
    end
end


for i=1:length(tbls)
    lst=find(~ismember(flds,tbls{i}.Properties.VariableNames));    
    for j=1:length(lst)
        if(ismember(types_safe{lst(j)},{"uint64","unit32","single","double"}))
            tmpdata=repmat(NaN,height(tbls{i}),1);
        elseif(types_safe{lst(j)}=="char")
            tmpdata=repmat(NaN,height(tbls{i}),1);
        else
            tmpdata=repmat({NaN},height(tbls{i}),1);
        end
        tbls{i}=[tbls{i} table(tmpdata,'VariableNames',{flds{lst(j)}})];
    end
end

for i=1:length(tbls)
    for j=1:length(types_safe)
        if(~ismember(types_safe{j},{'double','single','unit32','uint64'}))
            if(~iscell(tbls{i}.(flds{j})))
                if(ischar(tbls{i}.(flds{j})))
                    tbls{i}.(flds{j})=cellstr(tbls{i}.(flds{j}));
                else
                    tmp=cell(height(tbls{i}),1);
                    for k=1:height(tbls{i}.(flds{j}))
                        tmp{k}=tbls{i}.(flds{j})(k);
                    end
                    tbls{i}.(flds{j})=tmp;
                end
            end
        end
    end
end
out=vertcat(tbls{:});