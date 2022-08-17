function [tbl,ROIstats] = roiAverage( data, R, names ,splitrois,weighted)

if(~iscell(R)); R={R}; end;

if(nargin<3 || isempty(names))
    for i=1:length(R)
        names{i}=num2str(i);
    end
end
if ischar( names )
    names = {names};
end

if(nargin<4 || isempty(splitrois))
    % put all roi time courses in a single object or split
    splitrois=false;
end

if(nargin<5 || isempty(weighted))
    weighted=false;
end

if(length(data)>1)
    [tbl,ROIstats]=nirs.util.roiAverage(  data(1), R, names ,splitrois,weighted);
    
    if(isa(tbl,'table'))
        
        tbl=[table(repmat(cellstr(num2str(1)),height(tbl),1),'VariableNames',{'FileIdx'}) tbl];
    end
    for i=2:length(data)
        if(~isempty(ROIstats))
            [t,ROIstats(i,1)]=nirs.util.roiAverage(  data(i), R, names ,splitrois);
        else
            t=nirs.util.roiAverage(  data(i), R, names ,splitrois);
        end
        if(isa(tbl,'table'))
            t=[table(repmat(cellstr(num2str(i)),height(t),1),'VariableNames',{'FileIdx'}) t];
        end
        tbl=[tbl; t];
    end
    return
end





%[link, ilink] = sortrows(link, {'type','source', 'detector'});
% 
% allSrc=unique(link.source);
% allDet=unique(link.detector);
% 
% Rall={}; namesAll={};
% cnt=1;
% for idx=1:length(R)
%  
%     if(~ismember('weight',R{idx}.Properties.VariableNames))    
%         R{idx}.weight=ones(height(R{idx}),1);
%         R{idx}.weight=R{idx}.weight./sum(R{idx}.weight);
%     end
%     if(~ismember('Name',R{idx}.Properties.VariableNames) & ismember('name',R{idx}.Properties.VariableNames))    
%         R{idx}.Name=R{idx}.name;
%         R{idx}.name=[];
%     end
%     if(~ismember('Name',R{idx}.Properties.VariableNames))    
%         R{idx}.Name=repmat({names{idx}},height(R{idx}),1);
%     end
%        % Deal with any NaN's    
%     if(any(all([isnan(R{idx}.detector) isnan(R{idx}.source)],2)))
%         [s,d]=meshgrid(allSrc,allDet);
%         n=repmat(cellstr(R{idx}.Name{1}),size(s(:)));
%         w=repmat(R{idx}.weight(1),size(s(:)));
%         
%         R{idx}=[R{idx}; table(s(:),d(:),w(:),n,'VariableNames',{'source','detector','weight','Name'})];
%     end
%     
%     
%     
%     src=R{idx}.source(isnan(R{idx}.detector));
%     weight=R{idx}.weight(isnan(R{idx}.detector));
%     Name=R{idx}.Name(isnan(R{idx}.detector));
%     
%     if(length(src)>0)
%         R{idx}=[R{idx}; table(kron(src,ones(length(allDet),1)),...
%             kron(ones(length(src),1),allDet),kron(weight,ones(length(allDet),1)),...
%             reshape(repmat(Name',length(allDet),1),[],1),...
%             'VariableNames',{'source','detector','weight','Name'})];
%     end
%     
%     det=R{idx}.detector(isnan(R{idx}.source));
%     weight=R{idx}.weight(isnan(R{idx}.source));
%     Name=R{idx}.Name(isnan(R{idx}.source));
%     if(length(det)>0)
%       
%         
%         R{idx}=[R{idx}; table(kron(ones(length(det),1),allSrc),...
%             kron(ones(length(allSrc),1),det),kron(ones(length(allSrc),1),weight),...
%             reshape(repmat(Name',length(allSrc),1),[],1),...
%             'VariableNames',{'source','detector','weight','Name'})];
%     end
%     
%     R{idx}((isnan(R{idx}.source) | isnan(R{idx}.detector)),:)=[];
%     
% 
%     % Table contains the Names for the ROI.
%     names2=unique(R{idx}.Name);
%     sIdx=find(ismember(R{idx}.Properties.VariableNames,'source'));
%     dIdx=find(ismember(R{idx}.Properties.VariableNames,'detector'));
%     wIdx=find(ismember(R{idx}.Properties.VariableNames,'weight'));
%     
%     
%     lst=~ismember([R{idx}.source R{idx}.detector],[link.source link.detector],'rows');
%     R{idx}(lst,:)=[];
%     R{idx}=unique(R{idx});
%     
%     for idx2=1:length(names2)
%         lst=ismember(R{idx}.Name,names2{idx2});
%         Rall{cnt}=R{idx}(lst,[sIdx dIdx]);
%         namesAll{cnt}=names2{idx2};
%         if(~isempty(wIdx))
%             ContVect{cnt}=R{idx}.weight(lst);
%         end
%         cnt=cnt+1;
%     end
%    
% end
% R=Rall;
% names=namesAll;
% 
% if(length(data)>1)
%     if(isa(data,'nirs.core.Data'))
%         tbl=nirs.core.Data;
%         tbl(1)=[];
%         for idx=1:length(data)
%             tbl=[tbl; nirs.util.roiAverage(data(idx),R,names,splitrois)'];
%         end
%     else
%         tbl=[];
%         for idx=1:length(data)
%             thistbl=nirs.util.roiAverage(data(idx),R,names,splitrois);
%             description = {data(idx).description};
%             fileIdx=idx;
%             thistbl = [repmat(table(fileIdx,description),height(thistbl),1) thistbl];
%             tbl = [tbl; thistbl];
%         end
%     end
%     return
% end
% 
% 
% 
% 
% 
%     % The region definition is a table, parse it to the contrast vector
%     types=unique(link.type);
%     
%     RNew=cell(length(R)*length(types),1);
%     NamesNew=cell(length(R)*length(types),1);
%     cnt=1;
%      for idx2=1:length(R)
%     for idx=1:length(types)
%        
%             if(iscell(types))
%                 RNew{cnt}=ismember(link,[R{idx2} table(repmat({types{idx}},height(R{idx2}),1),...
%                     'VariableNames',{'type'})]);
%                 NamesNew{cnt}=[names{idx2} ':' types{idx}];
%             else
%                 RNew{cnt}=ismember(link,[R{idx2} table(repmat([types(idx)],height(R{idx2}),1),...
%                     'VariableNames',{'type'})]);
%                 NamesNew{cnt}=[names{idx2} ':' num2str(types(idx))];
%             end
%             
%             if(exist('ContVect'))
%                 ContVectNew{cnt}=ContVect{idx2};
%             end
%             cnt=cnt+1;
%         end
%     end
%     R=RNew;
%     namesOld=names;
%     names=NamesNew;
%     ContVect=ContVectNew;
    

if(isa(data,'nirs.core.Data'))
    data=data.sorted({'source','detector','type'});

    %First deal with the NaN values;
    [ContVect,R,names,namesOld,types]=nirs.util.ROIhelper(data,R,names);

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
            c(:,i)=c(:,i)/length(find(c(:,i)~=0));
        end
        if(weighted)
            W=diag(1./std(data.data,[],2));
            c=W*c;
            for i=1:size(c,2)
                 c(:,i)=c(:,i)/length(find(c(:,i)~=0));
            end
        end
        
        d(cnt)=nirs.core.Data;
        d(cnt).description=['ROI average' namesOld{floor((idx-1)/length(types))+1}];
        d(cnt).probe=nirs.core.Probe(NaN(1,3),NaN(1,3),table(repmat(1,length(types),1),...
            repmat(1,length(types),1),types,'VariableNames',{'source','detector','type'}));
       
        d(cnt).data=nan(size(data.data,1),size(c,2));
       for ii=1:size(d(cnt).data,1) 
            lst=find(~isnan(data.data(ii,:)));
            if(~isempty(lst))
                d(cnt).data(ii,:) = data.data(ii,lst)*c(lst,:);
            end
       end
       d(cnt).time=data.time;
       d(cnt).stimulus=data.stimulus;
       d(cnt).demographics=data.demographics;
       d(cnt).auxillary=data.auxillary;
       d(cnt).probe=nirs.core.ProbeROI({names{idx:idx+length(types)-1}});
       
       d(cnt).probe.fixeddistances=repmat(data.probe.distances(find(c(:,1)))/length(find(c(:,1))),size(c,2),1);
        cnt=cnt+1;
        
        
        
    end
    
    if(~splitrois)
        %combine the time courses into a single object
        dd=d(1);
      fixeddistances=d(1).probe.distances;
        for i=2:length(d)
            dd.data=[dd.data d(i).data];
            dd.probe.RegionNames={dd.probe.RegionNames{:} d(i).probe.RegionNames{:}}; 
            fixeddistances=[fixeddistances; d(i).probe.distances];
        end
        dd.probe.RegionNames=names;
        dd.probe.fixeddistances=fixeddistances;
        d=dd;
    end
    
    tbl=d;
    ROIstats=[];
    
elseif(isa(data,'nirs.core.sFCStats'))
     % loop over conditions
    data=data.sorted({'source', 'detector', 'type'}); 
     
    %First deal with the NaN values;
    [ContVect,R,names,namesOld,types]=nirs.util.ROIhelper(data,R,names);

    
    varnames = {'ROI_To','ROI_From', 'Contrast', 'R','Z','SE', 'DF', 'T', 'p'};
    tbl = table;

    Z = data.Z;
    
    
    vars=data.probe.link;
    %[vars, ivars] = sortrows(data.probe.link, {'source', 'detector', 'type'});
    %Z = data.Z(ivars,ivars,:);
    if(~isempty(data.ZstdErr))
        %covZ = data.ZstdErr(ivars, ivars,:).^2;
        covZ = data.ZstdErr.^2;
    else
        covZ=[];
    end
    
    
     % unique conditions
    uconds = unique(data.conditions, 'stable');

    
    % change ROIs to sorted indices
%     for i = 1:length(R)
%         R{i} = ilink(R{i});
%     end
    if(~exist('ContVect'))
        for i = 1:length(R)
            ContVect{i} = 1/length(find(R{i}));
        end
    end
    
    typ=zeros(length(names),1);
   
    
    
    for i = 1:length(uconds)
        
        b = Z(:,:,i);
        if(isempty(covZ))
            C=eye(size(Z,1));
        else
            C=covZ(:,:,i);
        end
        
        
        for j = 1:length(R)
            for j2=1:length(R)
                if(1) %j~=j2)
                    if(typ(j)==typ(j2))
                        
                        % contrast vector
                        c = zeros(size(b));
                        c(R{j},R{j2}) = ContVect{j}*ContVect{j2}';
                        %c=c/sum(c(:));
                        
                        CC=C;
                        CC(isnan(CC))=1E6;
%                          c2=c-diag(diag(c));  %remove the self terms
%                          c2= c2.*(triu(ones(size(c2))));
                        b(isnan(b))=1E6;
                        broi    = c(:)'*b(:)/sum(c(:));
                        rroi=tanh(broi);
                        df  =data.dfe;
                        
                        if(~isempty(covZ))
                            se = sqrt(c(:)'*diag(CC(:))*c(:));
                           % se = sqrt(1 - rroi.^2./(df-2));
                        else
                            se = sqrt((1 - rroi.^2)./(df-2));
                        end
                        
                        t       = rroi ./ se;
                        
                        
                        
                        
                        p       = 2*tcdf(-abs(t),df-2);
                        
                        tmp = cell2table({names{j},names{j2} uconds{i}, rroi,broi, se, df, t, p});
                        tmp.Properties.VariableNames = varnames;
                        
                        
                        tbl = [tbl; tmp];
                    end
                end
            end
        end
    end
    q   = nirs.math.fdr( tbl.p );
    [~,power] = nirs.math.MDC(tbl,.8,.05);
    tbl = [tbl table(q) table(power)];
    
    ROIstats=[];
    if(nargout==2)
        warning('havent created roi connectivity model yet');
    end
    tbl(isnan(tbl.power),:)=[];
    
elseif(isa(data,'nirs.core.ChannelFStats'))
    data=data.sorted({'cond', 'source', 'detector', 'type'});
    
        %First deal with the NaN values;
    [ContVect,R,names,namesOld,types]=nirs.util.ROIhelper(data,R,names);

%    [vars, ivars] = sortrows(data.variables, {'cond', 'source', 'detector', 'type'});
    vars=data.variables;
 
    F = data.F; %(ivars);
    df1 = data.df1; %(ivars);
    df2 = data.df2; %(ivars);
    uconds = unique(vars.cond, 'stable');
   
  
   % loop over conditions
    varnames = {'ROI','type', 'Contrast','F', 'DF1', 'DF2', 'p'};
    tbl = table;
    
     % change ROIs to sorted indices
%     for i = 1:length(R)
%         R{i} = ilink(R{i});
%     end
    if(~exist('ContVect'))
        for i = 1:length(R)
            ContVect{i} = 1/length(find(R{i}));
        end
    end
    
    cc=zeros(size(F,1),length(R)*length(uconds));
    vvs =table;
    for i = 1:length(uconds)
        lst = strcmp(vars.cond, uconds{i});
        b = F(lst);
        d1 = df1(lst);
        d2 = df2(lst);
        
        for j = 1:length(R)
            % contrast vector
            c = zeros(size(b));
            c(R{j}) = ContVect{j};
           
            
            cc(lst,(i-1)*length(R)+j)=c;
            vvs = [vvs; table(namesOld(floor((j-1)/length(types))+1), types(mod(j-1,length(types))+1),uconds(i),'VariableNames',{'ROI','type','cond'})];
            froi    = c'*b;
            df1roi     = c'*d1;
            df2roi     = c'*d2;
            p       =  fcdf( 1./froi, df2roi,df1roi);
            %p       =  fcdf( 1./froi, df1roi,df2roi);
            
            tmp = cell2table({namesOld(floor((j-1)/length(types))+1),...
                types(mod(j-1,length(types))+1), uconds{i},  froi, df1roi,df2roi,p});
            tmp.Properties.VariableNames = varnames;
            
            
            tbl = [tbl; tmp];
        end
    end
    
    q   = nirs.math.fdr( tbl.p );
    tbl = [tbl table(q)];
    
    ROIstats=nirs.core.ChannelFStats;
    ROIstats.description='region of interest F-stats';
    ROIstats.F=tbl.F;
    ROIstats.df1=tbl.DF1;
    ROIstats.df2=tbl.DF2;
    ROIstats.demographics=data.demographics;
    ROIstats.variables=vvs;
    
    ROIstats.probe=nirs.core.ProbeROI(names);
    
    
    
else
    
    % sort variables
    data=data.sorted({'cond', 'source', 'detector', 'type'});
    
        %First deal with the NaN values;
    [ContVect,R,names,namesOld,types]=nirs.util.ROIhelper(data,R,names);

    
    beta = data.beta; %(ivars);
    covb = data.covb; %(ivars, ivars);
    vars=data.table;
     % unique conditions
    
    uconds = unique(vars.cond, 'stable');
    
    % loop over conditions
    varnames = {'ROI','type', 'Contrast','Beta', 'SE', 'DF', 'T', 'p'};
    tbl = table;
    
    
    % change ROIs to sorted indices
%     for i = 1:length(R)
%         R{i} = ilink(R{i});
%     end
    if(~exist('ContVect'))
        for i = 1:length(R)
            ContVect{i} = 1/length(find(R{i}));
        end
    end
    
    cc=zeros(size(beta,1),length(R)*length(uconds));
    vvs =table;
    models={};
    for j = 1:length(R)
        for i = 1:length(uconds)
            lst = strcmp(vars.cond, uconds{i});
            b = beta(lst);
            C = covb(lst,lst);
            
            
            
            % contrast vector
            c = zeros(size(b));
            c(R{j}) = ContVect{j};
            
            if(weighted)
                W=diag(1./diag(data.covb(lst,lst)));
                
                c=W*c;
                c=c./sum(c);
            end
            
            broi=nan(size(c,2),size(b,2));
            se=nan(size(c,2),size(c,2));
            
            cc(lst,(i-1)*length(R)+j)=c;
            vvs = [vvs; table(namesOld(floor((j-1)/length(types))+1), types(mod(j-1,length(types))+1),uconds(i),'VariableNames',{'ROI','type','cond'})];
            lst=find(~isnan(b));
            %             broi(lst)    = c(lst,:)'*b(lst,:);
            %             se(lst)      = sqrt(c(lst,:)'*C(lst,lst)*c(lst,:));
            
            c=c*sum(c)/sum(c(lst)); %---- Ted: need to check here for generalizability
            
            broi    = c(lst,:)'*b(lst,:); %---- Ted: check here
            se      = diag(sqrt(c(lst,:)'*C(lst,lst)*c(lst,:))); %---- Ted: check here
            
            t       = broi./se;
            df      = data.dfe;
            p       = 2*tcdf(-abs(t),df);
            
            tmp = cell2table({namesOld(floor((j-1)/length(types))+1),...
                types(mod(j-1,length(types))+1), uconds{i},  broi, se, df, t, p});
            tmp.Properties.VariableNames = varnames;
            
            tbl = [tbl; tmp];
        end
        
        if(ismember('model',vars.Properties.VariableNames))
            % include the LinearModel (diagnotics) in the ROI
            str=matlab.lang.makeValidName(uconds);
            % 
            model = combineLinearModels(repmat(c,length(uconds),1),vars.model,data.WhiteningW,str);
            tmp.model={model};
            t=model.Coefficients;
            
            [~,l]=ismember(t.Properties.RowNames,str);
               
            lst2=find(ismember(tbl.ROI,tmp.ROI) & ismember(tbl.type,tmp.type));
            
            % Use the values from the table instead
%             tbl(lst2,:).Beta=t.Estimate(l);
%             tbl(lst2,:).SE=t.SE(l);
%             tbl(lst2,:).T=t.tStat(l);
%             tbl(lst2,:).p=t.pValue(l);
%             tbl(lst2,:).DF=model.DFE*ones(length(l),1);
            for i=1:length(uconds)
                models{end+1,1}=model;
            end
        end
        
    end
    
    q   = nirs.math.fdr( tbl.p );
    [~,power] = nirs.math.MDC(tbl,.8,.05);
    tbl = [tbl table(q) table(power)];
    
    
    ROIstats=nirs.core.ChannelStats;
    ROIstats.description='region of interest stats';
    ROIstats.beta=tbl.Beta;
    ROIstats.dfe=tbl.DF;
    ROIstats.covb=cc'*covb*cc; %----- Ted: need to check this line for compatability with NaNs
    ROIstats.demographics=data.demographics;
    ROIstats.basis=data.basis;
    
    if(ismember('model',vars.Properties.VariableNames))
        tbl=[tbl table(models,'VariableName',{'model'})];
    end
   
    
    if(ismember('model',vars.Properties.VariableNames))
        ROIstats.categoricalvariableInfo=data.categoricalvariableInfo;
        vvs=[vvs table(tbl.model,'VariableName',{'model'})];
    end
    
    ROIstats.variables=vvs;
    
    ROIstats.probe=nirs.core.ProbeROI(names);
    tbl(isnan(tbl.power),:)=[];    
end



end



function mdl = combineLinearModels(c,models,WhiteningW,str)

tbl=table();
w=[];

if(~isempty(WhiteningW))
     iW=inv(WhiteningW);
else
   iW=1;
end

tbl=models{1}.Variables;
flds=tbl.Properties.VariableNames;
for j=1:length(flds)
    tbl.(flds{j})(:)=0;
end

for idx=1:length(models);
    tmp=models{idx}.Variables;
    for j=1:length(flds)
        tbl.(flds{j})=tbl.(flds{j})+tmp.(flds{j})*c(idx);
    end
end

flds=tbl.Properties.VariableNames;
mask=zeros(height(tbl),1);
for i=1:length(flds)
    if(~strcmp(flds{i},models{1}.ResponseName))
        mask=max(mask,abs(tbl.(flds{i})));
    end
end

lst=find(mask<1E-9);
tbl(lst,:)=[];

mdl=fitlm(tbl,models{1}.Formula,'dummyVarCoding','full');

%     
%     if(c(idx)~=0)
%         thistbl=models{idx}.Variables;
%         flds=thistbl.Properties.VariableNames;
%         for j=1:length(flds)
%             if(~ismember(flds{j},str))
%                 thistbl.(flds{j})=iW*thistbl.(flds{j})*c(idx);
%             else
%                 thistbl.(flds{j})=iW*thistbl.(flds{j});
%             end
%         end
%        %thistbl.beta=thistbl.beta*c(idx);
% 
%         tbl=[tbl; thistbl];
%         w=[w;  models{idx}.ObservationInfo.Weights];
%     end
% end
% 
% flds=tbl.Properties.VariableNames;
% mask=zeros(height(tbl),1);
% for i=1:length(flds)
%     if(~strcmp(flds{i},models{1}.ResponseName))
%         mask=max(mask,abs(tbl.(flds{i})));
%     end
% end
% 
% lst=find(mask<1E-9);
% tbl(lst,:)=[];
% w(lst)=[];
% 
% 
% mdl=fitlm(tbl,models{1}.Formula,'weights',max(w,eps(1)),'dummyVarCoding','full');

end
