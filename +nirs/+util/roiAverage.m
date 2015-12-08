function tbl = roiAverage( data, R, names )



if(~iscell(R) && isa(R,'table'))
    if(any(ismember(R.Properties.VariableNames,'Name')))
        % Table contains the Names for the ROI.
        names=unique(R.Name);
        sIdx=find(ismember(R.Properties.VariableNames,'source'));
        dIdx=find(ismember(R.Properties.VariableNames,'detector'));
        wIdx=find(ismember(R.Properties.VariableNames,'weight'));
        
        for idx=1:length(names)
            lst=ismember(R.Name,names{idx});
            RR{idx}=R(lst,[sIdx dIdx]);
            if(~isempty(wIdx))
                ContVect{idx}=R.weight(lst);   
            end
            
        end
        R=RR;
    end
end
   
    if(nargin<3 && ~exist('names'))
        names=cellstr(num2str(1:length(R)));
    end
    if ischar( names )
        names = {names};
    end 
    
    if(length(data)>1)
        if(isa(data,'nirs.core.Data'))
            tbl=nirs.core.Data;
            tbl(1)=[];
            for idx=1:length(data)
                tbl=[tbl; nirs.util.roiAverage(data(idx),R,names)'];
            end
        else
            tbl=[];
            for idx=1:length(data)
                thistbl=nirs.util.roiAverage(data(idx),R,names);
                description = {data(idx).description};
                fileIdx=idx;
                thistbl = [repmat(table(fileIdx,description),height(thistbl),1) thistbl];
                tbl = [tbl; thistbl];
            end
        end
        return
    end
    
    data=data.sorted;
     % sort probe
    link = data.probe.link;
    [link, ilink] = sortrows(link, {'type','source', 'detector'});
    
    if(isa(R{1},'table'))
        
        %First deal with the NaN values;
        allSrc=unique(link.source);
        allDet=unique(link.detector);
        
        for idx=1:length(R)
            if(any(all([isnan(R{idx}.detector) isnan(R{idx}.source)],2)))
                [s,d]=meshgrid(allSrc,allDet);
                R{idx}=[R{idx}; table(s(:),d(:),'VariableNames',{'source','detector'})];
            end
            
            
            src=R{idx}.source(isnan(R{idx}.detector));
            if(length(src)>0)
            R{idx}=[R{idx}; table(kron(src,ones(length(allDet),1)),...
                kron(ones(length(src),1),allDet),'VariableNames',{'source','detector'})];
            end
            
            det=R{idx}.detector(isnan(R{idx}.source));
            if(length(det)>0)
            R{idx}=[R{idx}; table(kron(ones(length(det),1),allSrc),...
                    kron(ones(length(allSrc),1),det),'VariableNames',{'source','detector'})];
            end
            
        end
        
        
        % The region definition is a table, parse it to the contrast vector
        types=unique(link.type);
        
        RNew=cell(length(R)*length(types),1);
        NamesNew=cell(length(R)*length(types),1);
        cnt=1;
        for idx=1:length(types)
            for idx2=1:length(R)
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
        if(exist('ContVect'))
            ContVect=ContVectNew;
        end
        
    end
    
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
                end
                
                tbl = [tbl; tmp];
            end
        end
        
        q   = nirs.math.fdr( tbl.p );
        tbl = [tbl table(q)];
    end
end


function mdl = combineLinearModels(c,models)

tbl=table();
w=[];

for idx=1:length(models);
    if(c(idx)~=0)
        tbl=[tbl; models{idx}.Variables];
        w=[w; c(idx)*models{idx}.ObservationInfo.Weights];
    end
end
mdl=fitlm(tbl,models{1}.Formula,'weights',w);


end
