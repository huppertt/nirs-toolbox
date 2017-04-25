function tbl = roiAverage( data, R, names )

if(~iscell(R)); R={R}; end;

if(nargin<3 && ~exist('names'))
    for i=1:length(R)
        names{i}=num2str(i);
    end
end
if ischar( names )
    names = {names};
end


%First deal with the NaN values;
data=data.sorted;
% sort probe
link = data(1).probe.link;
[link, ilink] = sortrows(link, {'type','source', 'detector'});

allSrc=unique(link.source);
allDet=unique(link.detector);

Rall={}; namesAll={};
cnt=1;
for idx=1:length(R)
 
    if(~ismember('weight',R{idx}.Properties.VariableNames))    
        R{idx}.weight=ones(height(R{idx}),1);
    end
    if(~ismember('Name',R{idx}.Properties.VariableNames) & ismember('name',R{idx}.Properties.VariableNames))    
        R{idx}.Name=R{idx}.name;
        R{idx}.name=[];
    end
    if(~ismember('Name',R{idx}.Properties.VariableNames))    
        R{idx}.Name=repmat({names{idx}},height(R{idx}),1);
    end
       % Deal with any NaN's    
    if(any(all([isnan(R{idx}.detector) isnan(R{idx}.source)],2)))
        [s,d]=meshgrid(allSrc,allDet);
        n=repmat(cellstr(R{idx}.Name{1}),size(s(:)));
        w=repmat(R{idx}.weight(1),size(s(:)));
        
        R{idx}=[R{idx}; table(s(:),d(:),w(:),n,'VariableNames',{'source','detector','weight','Name'})];
    end
    
    
    
    src=R{idx}.source(isnan(R{idx}.detector));
    weight=R{idx}.weight(isnan(R{idx}.detector));
    Name=R{idx}.Name(isnan(R{idx}.detector));
    
    if(length(src)>0)
        R{idx}=[R{idx}; table(kron(src,ones(length(allDet),1)),...
            kron(ones(length(src),1),allDet),kron(weight,ones(length(allDet),1)),...
            reshape(repmat(Name',length(allDet),1),[],1),...
            'VariableNames',{'source','detector','weight','Name'})];
    end
    
    det=R{idx}.detector(isnan(R{idx}.source));
    weight=R{idx}.weight(isnan(R{idx}.source));
    Name=R{idx}.Name(isnan(R{idx}.source));
    if(length(det)>0)
      
        
        R{idx}=[R{idx}; table(kron(ones(length(det),1),allSrc),...
            kron(ones(length(allSrc),1),det),kron(ones(length(allSrc),1),weight),...
            reshape(repmat(Name',length(allSrc),1),[],1),...
            'VariableNames',{'source','detector','weight','Name'})];
    end
    
    R{idx}((isnan(R{idx}.source) | isnan(R{idx}.detector)),:)=[];
    

    % Table contains the Names for the ROI.
    names2=unique(R{idx}.Name);
    sIdx=find(ismember(R{idx}.Properties.VariableNames,'source'));
    dIdx=find(ismember(R{idx}.Properties.VariableNames,'detector'));
    wIdx=find(ismember(R{idx}.Properties.VariableNames,'weight'));
    
    
    lst=~ismember([R{idx}.source R{idx}.detector],[link.source link.detector],'rows');
    R{idx}(lst,:)=[];
    R{idx}=unique(R{idx});
    
    for idx2=1:length(names2)
        lst=ismember(R{idx}.Name,names2{idx2});
        Rall{cnt}=R{idx}(lst,[sIdx dIdx]);
        namesAll{cnt}=names2{idx2};
        if(~isempty(wIdx))
            ContVect{cnt}=R{idx}.weight(lst);
        end
        cnt=cnt+1;
    end
   
end
R=Rall;
names=namesAll;

if(length(data)>1)
    if(isa(data,'nirs.core.Data'))
        tbl=nirs.core.Data;
        tbl(1)=[];
        for idx=1:length(data)
            tbl=[tbl; nirs.util.roiAverage(data(idx),R)'];
        end
    else
        tbl=[];
        for idx=1:length(data)
            thistbl=nirs.util.roiAverage(data(idx),R);
            description = {data(idx).description};
            fileIdx=idx;
            thistbl = [repmat(table(fileIdx,description),height(thistbl),1) thistbl];
            tbl = [tbl; thistbl];
        end
    end
    return
end





    % The region definition is a table, parse it to the contrast vector
    types=unique(link.type);
    
    RNew=cell(length(R)*length(types),1);
    NamesNew=cell(length(R)*length(types),1);
    cnt=1;
     for idx2=1:length(R)
    for idx=1:length(types)
       
            if(iscell(types))
                RNew{cnt}=ismember(link,[R{idx2} table(repmat({types{idx}},height(R{idx2}),1),...
                    'VariableNames',{'type'})]);
                NamesNew{cnt}=[names{idx2} ':' types{idx}];
            else
                RNew{cnt}=ismember(link,[R{idx2} table(repmat([types(idx)],height(R{idx2}),1),...
                    'VariableNames',{'type'})]);
                NamesNew{cnt}=[names{idx2} ':' num2str(types(idx))];
            end
            
            if(exist('ContVect'))
                ContVectNew{cnt}=ContVect{idx2};
            end
            cnt=cnt+1;
        end
    end
    R=RNew;
    namesOld=names;
    names=NamesNew;
    ContVect=ContVectNew;
    

if(isa(data,'nirs.core.Data'))
    data=data.sorted({'type','source','detector'});
    if(~exist('ContVect'))
        for i = 1:length(R)
            ContVect{i} = 1/length(find(R{i}));
        end
    end
    clear d; cnt=1;
    for idx=1:length(types):length(R)
        
        c = zeros(size(data.data,2),length(types));
        for i=1:length(types)
            c(R{idx+i-1},i) = ContVect{idx+i-1};
            c(:,i)=c(:,i)/sum(c(:,i));
        end
        d(cnt)=nirs.core.Data;
        d(cnt).description=['ROI average' namesOld{floor(idx/length(types))+1}];
        d(cnt).probe=nirs.core.Probe(NaN(1,3),NaN(1,3),table(repmat(1,length(types),1),...
            repmat(1,length(types),1),types,'VariableNames',{'source','detector','type'}));
        d(cnt).data = data.data*c;
        d(cnt).time=data.time;
        d(cnt).stimulus=data.stimulus;
        cnt=cnt+1;
    end
    tbl=d;
    
else
    
    % sort variables
    [vars, ivars] = sortrows(data.variables, {'cond', 'source', 'detector', 'type'});
    beta = data.beta(ivars);
    covb = data.covb(ivars, ivars);
    
    % unique conditions
    uconds = unique(vars.cond, 'stable');
    
    % loop over conditions
    varnames = {'ROI', 'Contrast', 'Beta', 'SE', 'DF', 'T', 'p'};
    tbl = table;
    
    
    % change ROIs to sorted indices
    for i = 1:length(R)
        R{i} = ilink(R{i});
    end
    if(~exist('ContVect'))
        for i = 1:length(R)
            ContVect{i} = 1/length(find(R{i}));
        end
    end
    
    for i = 1:length(uconds)
        lst = strcmp(vars.cond, uconds{i});
        b = beta(lst);
        C = covb(lst,lst);
        
        for j = 1:length(R)
            % contrast vector
            c = zeros(size(b));
            c(R{j}) = ContVect{j};
            c=c/sum(c);
           
            broi    = c'*b;
            se      = sqrt(c'*C*c);
            t       = broi / se;
            df      = data.dfe;
            p       = 2*tcdf(-abs(t),df);
            
            
            tmp = cell2table({names{j}, uconds{i}, broi, se, df, t, p});
            tmp.Properties.VariableNames = varnames;
            
            if(ismember('model',vars.Properties.VariableNames))
                % include the LinearModel (diagnotics) in the ROI
                model = combineLinearModels(c,vars.model(lst));
                tmp.model={model};
                t=model.Coefficients;
                str=model.PredictorNames{1};
                l=find(ismember(t.Properties.RowNames,[str '_' tmp.Contrast{1}]));
                
                % Use the values from the table instead
                tmp.Beta=t.Estimate(l);
                tmp.SE=t.SE(l);
                tmp.T=t.tStat(l);
                tmp.p=t.pValue(l);
                tmp.DF=model.DFE;
                
            end
            
            tbl = [tbl; tmp];
        end
    end
    
    q   = nirs.math.fdr( tbl.p );
    [~,power] = nirs.math.MDC(tbl,.8,.05);
    tbl = [tbl table(q) table(power)];
end


end



function mdl = combineLinearModels(c,models)

tbl=table();
w=[];
c=c/rms(c);

for idx=1:length(models);
    if(c(idx)~=0)
        thistbl=models{idx}.Variables;
        thistbl.beta=thistbl.beta*c(idx);
        tbl=[tbl; thistbl];
        w=[w; models{idx}.ObservationInfo.Weights];
    end
end
mdl=fitlm(tbl,models{1}.Formula,'weights',max(w,eps(1)),'dummyVarCoding','full');


end
