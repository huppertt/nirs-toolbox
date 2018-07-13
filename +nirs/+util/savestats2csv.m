function varargout = savestats2csv(S,file,weight)
% This function saves all the Subject level stats to a CSV file for import
% into SPSS/SAS

% % Make sure the data is already in subject-level
% j=nirs.modules.SubjLevelStats
% j.sortfield='subject';
% S=j.run(SubjStats);

if(nargin<3)
    weight='none';
end


tblall={};

% Now create the data table
for i=1:length(S)
    demo=nirs.createDemographicsTable(S(i));
        
    for j=1:length(demo.Properties.VariableNames)
        if(isstr(demo.(demo.Properties.VariableNames{j})))
            demo.(demo.Properties.VariableNames{j})=...
                cellstr(demo.(demo.Properties.VariableNames{j}));
        end
    end
    
    C = S(i).covb;
    [u,s,v]=nirs.math.mysvd(C);
    W = inv(chol(C));
    
    tbl=S(i).table;
    tbl=sortrows(tbl,{'type','cond','source','detector'});
    
    columns = strcat(repmat(cellstr('Src'),height(tbl),1),num2str(tbl.source),...
        repmat(cellstr('_Det'),height(tbl),1),num2str(tbl.detector),...
        repmat(cellstr('_'),height(tbl),1),tbl.cond,repmat(cellstr('_'),height(tbl),1),...
        tbl.type);
    
    b=tbl.beta;
    if(strcmp(weight,'whiten'))
        b = W*b;
    elseif('stderr')
        b=[b; 1./tbl.se];
        cw=strcat(columns,repmat(cellstr('_weight'),length(columns),1));
        columns={columns{:} cw{:}}';
    end
    BB={};
    for j=1:length(b); BB{j}=b(j); end;
    tbl=cell2table(BB);
    tbl.Properties.VariableNames=columns;
    
    tblall{i} =[demo tbl];
   % tblall=[tblall; tbllocal];
end
    
n={};
for i=1:length(tblall);
    n={n{:} tblall{i}.Properties.VariableNames{:}};
end
n=unique(n);
for i=1:length(tblall)
    lst=find(~ismember(n,tblall{i}.Properties.VariableNames));
    tblall{i}=[tblall{i} cell2table(cell(1,length(lst)),'VariableNames',{n{lst}})];
end

tblall=vertcat(tblall{:});


if(nargout==1)
    varargout{1}=tblall;
end

if(nargin>1 & ~isempty(file))
    writetable(tblall,file);
end
    
    
   


