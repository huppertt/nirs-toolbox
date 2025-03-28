function f=draw( obj, vtype, vrange, thresh ,flip)
%% draw - Draws channelwise values on a probe.
% Args:
%     vtype   - either 'Pearsons' or 'Fischer-Z' or 'Grangers' or 'Grangers-F'
%     vrange  - range to display; either a scalar or vector with 2 elements
%     thresh  - either a scalar such that values > thresh are significant or
%               a string specifying statistical significance (e.g. 'p < 0.05')
%
% Examples:
%     stats.draw( 'Pearsons', [-5 5], 'p < 0.1' )
%     stats.draw( 'Grangers', 5, 3 )

% type is either beta or tstat
if nargin < 2;
    vtype='R';
    
    %         if(strcmp(obj.type,'Grangers'));
    %              vtype = 'Grangers';
    %         else
    %             vtype='R';
    %
    %         end
end
if(~isempty(strfind(vtype,'matrix')))
    vtype(strfind(vtype,'matrix')-1+[1:length('matrix')])=[];
    drawtype='matrix';
else
    drawtype='line';
end

if(~isempty(strfind(vtype,'line')))
    vtype(strfind(vtype,'line')-1+[1:length('line')])=[];
end

if nargin < 3
    vrange_arg = [];
else
    vrange_arg = vrange;
end

vtype(strfind(vtype,' '))=[];

if(nargin<5 || isempty(flip))
    flip=[1 1];
end




for cIdx=1:length(obj.conditions)
    tbl=obj.table;
    tbl(~ismember(tbl.condition,obj.conditions{cIdx}),:)=[];
    
    switch(lower(vtype))
        %         case('grangers')
        %             values=tbl.Grangers;
        %             pval = tbl.pvalue;
        %             cmap=flipdim(colormap(autumn(2001)),1);
        %         case('grangers-f')
        %             values=tbl.F;
        %             pval = tbl.pvalue;
        %              cmap=flipdim(colormap(autumn(2001)),1);
        case('r')
            values=tbl.R;
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
            clabel = 'r-value';
        case('z')
            values=tbl.Z;
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
            clabel = 'Z';
        case('t')
            values=tbl.t;
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
            clabel = 't-statistic';
        otherwise
            warning('type not recognized: using Correlation');
            values=tbl.R;
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
            %cmap=flipdim(colormap(autumn(2001)),1);
    end

    pval=tbl.pvalue;
    qval=tbl.qvalue;
    

    if(isa(obj.probe,'nirs.core.ProbeHyperscan'))

        % mask within subject
        mask=ones(height(obj.probe.link));
        for i=1:length(obj.probe.SubjectLabels)
            lstt=find(ismember(obj.probe.link.SubjectLabel,obj.probe.SubjectLabels{i}));
            mask(lstt,lstt)=0;
        end
        values(find(mask==0))=NaN;
        pval(find(mask==0))=NaN;
        qval(find(mask==0))=NaN;
        %qval=pval;
        %lstt=find(~isnan(pval));
        %qval(lstt)=nirs.math.BenjaminiHochberg(pval(lstt));
    end

    
    % significance mask
    if nargin < 4 || isempty(thresh)
        mask = ~isnan(values);
        
    elseif isscalar(thresh)
        mask = abs(values) > thresh;
        
    elseif isvector(thresh) && isnumeric(thresh)
        mask = values < thresh(1) | values > thresh(2);
        
    elseif isstr(thresh)
        % takes in string in form of 'p < 0.05' or 'q < 0.10'
        s = strtrim( strsplit( thresh, '<' ) );
        
        if(~isempty(strfind(s{1},'q')))
            mask = qval < str2double(s{2});
        else
            mask = pval < str2double(s{2});
        end
    end
    I=eye(sqrt(length(mask)));
    mask=mask.*(~I(:));
    
    % range to show
    if isempty(vrange_arg)
        vmax    = max(abs(values(:).*mask(:)));
        vrange  = vmax*[-1 1];
    else
        vrange  = vrange_arg;
    end
    
    p = obj.probe;
    
    % Expand ROI link to channel
    if iscell(p.link.source)
        link = p.link;
        link2 = table([],[],{},'VariableNames',{'source','detector','type'});
        hyper = '';
        for i=1:height(link)
            for j=1:length(link.source{i})
                link2(end+1,:) = table(link.source{i}(j),link.detector{i}(j),link.type(i));
                if ismember('hyperscan',link.Properties.VariableNames)
                    hyper(end+1,1) = link.hyperscan(i);
                end
            end
        end
        if ismember('hyperscan',link.Properties.VariableNames)
            link2.hyperscan = hyper;
        end
        p.link = link2;
    end
    
    typesOrigin=tbl.TypeOrigin;
    typesDest=tbl.TypeDest;
    
    
    % convert to strings for consistency in loop below
    if any(isnumeric( typesOrigin))
        typesOrigin = cellfun(@(x) {num2str(x)}, num2cell( typesOrigin));
    end
    
    % convert to strings for consistency in loop below
    if any(isnumeric(typesDest))
        typesDest = cellfun(@(x) {num2str(x)}, num2cell(typesDest));
    end
    
    
    % unique types
    utypesOrigin = unique(typesOrigin, 'stable');
    utypesDest = unique(typesDest, 'stable');
    
    
    % colormap
    
    z = linspace(vrange(1), vrange(2), size(cmap,1))';
    
    f(cIdx)=figure;
    
   
        
    
        cnt=1;
        for ii=1:length(utypesOrigin)
            for jj=1:length(utypesDest)
                if(strcmp(utypesOrigin(ii),utypesDest(jj)))
                    
                    lst=find(strcmp(typesOrigin,utypesOrigin(ii)) & ...
                        strcmp(typesDest,utypesDest(jj)));
                    
                    
                    vals = values(lst);
                    
                    % this mask
                    m = mask(lst);
                    
                    % map to colors
                    idx = bsxfun(@minus, vals', z);
                    [~, idx] = min(abs(idx), [], 1);
                    colors = cmap(idx, :);
                    
                    if(strcmp(drawtype,'line'))
                        figure(f(cIdx));
                        sp=subplot(length(utypesOrigin),1,cnt);
                        
                        h=obj.probe.draw([],[],sp);
                        set(h,'Color', [.7 .7 .7]);
                        axes(sp);
                        if(isa(obj.probe,'nirs.core.ProbeHyperscan1020') |...
                                isa(obj.probe,'nirs.core.ProbeHyperscan'))
                            srcPos=obj.probe.srcPos_drawing;
                            detPos=obj.probe.detPos_drawing;

                        else

                            link=obj.probe.link;
                            link=link(ismember(link.type,link.type(1)),:);
                            for id=1:length(h)
                                XYZ=[get(h(id),'XData')' get(h(id),'YData')' get(h(id),'ZData')'];
                                if(isempty(get(h(id),'ZData')))
                                    XYZ(:,3)=0;
                                end
                                srcPos(link.source(id),:)=XYZ(1,:);
                                detPos(link.detector(id),:)=XYZ(2,:);

                            end
                        end
                        
                        posOrig=(srcPos(tbl.SourceOrigin(lst),:)+...
                            detPos(tbl.DetectorOrigin(lst),:))/2;
                        posDest=(srcPos(tbl.SourceDest(lst),:)+...
                            detPos(tbl.DetectorDest(lst),:))/2;
                        
                        
                        X=[posOrig(:,1) posDest(:,1)];
                        Y=[posOrig(:,2) posDest(:,2)];
                        Z=[posOrig(:,3) posDest(:,3)];
                        
                       
                        h2=[];
                        for idx=1:length(vals)
                            if(m(idx))
                                h2(end+1)=line(sp,X(idx,:),Y(idx,:),Z(idx,:),'Color',colors(idx,:),'linewidth',3);
                                set(h2(end),'Tag','fNIRS_ConnLine');
                            end
                        end
                        
                        colormap(cmap);
                        colorbar(sp);
                        caxis(vrange);
                       % title(sp,[utypesOrigin{ii} ' --> ' utypesDest{jj}], 'Interpreter','none')
                        
                    else
                        
                        figure(f(cIdx));
                        sp=subplot(length(utypesOrigin),1,cnt);
                        
                        LabelsOrig=strcat(repmat('src-',length(lst),1),num2str(tbl.SourceOrigin(lst)),...
                            repmat(':det-',length(lst),1), num2str(tbl.DetectorOrigin(lst)));
                        
                        LabelsDet=strcat(repmat('src-',length(lst),1),num2str(tbl.SourceDest(lst)),...
                            repmat(':det-',length(lst),1), num2str(tbl.DetectorDest(lst)));
                        [LabelsDet,i]=unique(LabelsDet,'rows');
                        [LabelsOrig,j]=unique(LabelsOrig,'rows');
                                                
                        vals=reshape(vals(:).*m,length(i),length(j));
                        imagesc(sp,vals,[vrange]);
                        colorbar;
                        set(sp,'XTick',[1:length(i)]);
                        set(sp,'YTick',[1:length(j)]);
                        
                        set(sp,'YTickLabel',{LabelsDet})
                        set(sp,'XTickLabel',{LabelsOrig})
                        set(sp,'XtickLabelRotation',90);
                        title(sp,[utypesOrigin{ii} ' --> ' utypesDest{jj}], 'Interpreter','none')
                    end
                    cnt=cnt+1;
                end
            end
            
        end
    
    set(f(cIdx),'Name',obj.conditions{cIdx},'NumberTitle','off')
    %supertitle(f(cIdx),obj.conditions{cIdx});
end

