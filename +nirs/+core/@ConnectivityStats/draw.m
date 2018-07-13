function draw( obj, vtype, vrange, thresh )
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
        if(strcmp(obj.type,'Grangers'));
             vtype = 'Grangers'; 
        else
            vtype='pearsons';
            
        end
    end
    
    tbl=obj.table;
    
    switch(lower(vtype))
        case('grangers')
            values=tbl.Grangers;
            pval = tbl.pvalue;
            cmap=flipdim(colormap(autumn(2001)),1);
        case('grangers-f')
            values=tbl.F;
            pval = tbl.pvalue;
             cmap=flipdim(colormap(autumn(2001)),1);
        case('pearsons')
            values=tbl.Pearsons;
            pval=tbl.pvalue;
             [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
        case('fischer-z')
            values=tbl.Z;
            pval=tbl.pvalue;
             [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
        otherwise
            warning('type not recognized: using Grangers');
             values=tbl.Grangers;
            pval = tbl.pvalue;
             cmap=flipdim(colormap(autumn(2001)),1);
    end
    
    
    % range to show
    if nargin < 3 || isempty(vrange)
        vmax    = max(abs(values(:)));
        vrange  = vmax*[-1 1];
    end

    % significance mask
    if nargin < 4
        mask = ones(size(values)) > 0;
        
    elseif isscalar(thresh)
        mask = abs(values) > thresh;
        
    elseif isvector(thresh) && isnumeric(thresh)
        mask = values < thresh(1) | values > thresh(2);
        
    elseif isstr(thresh)
        % takes in string in form of 'p < 0.05' or 'q < 0.10'
        s = strtrim( strsplit( thresh, '<' ) );
        
        mask = pval < str2double(s{2});
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
  
    f(1)=figure;
    f(2)=figure;
    
    cnt=1;
    for ii=1:length(utypesOrigin)
        for jj=1:length(utypesDest)
            lst=find(strcmp(typesOrigin,utypesOrigin(ii)) & ...
                strcmp(typesDest,utypesDest(jj)));
            
            
            vals = values(lst);
            
            % this mask
            m = mask(lst);
            
            % map to colors
            idx = bsxfun(@minus, vals', z);
            [~, idx] = min(abs(idx), [], 1);
            colors = cmap(idx, :);
            
            posOrig=(obj.probe.srcPos(tbl.SourceOrigin(lst),:)+...
                obj.probe.detPos(tbl.DetectorOrigin(lst),:))/2;
            posDest=(obj.probe.srcPos(tbl.SourceDest(lst),:)+...
                obj.probe.detPos(tbl.DetectorDest(lst),:))/2;
            
          
            X=[posOrig(:,1) posDest(:,1)];
            Y=[posOrig(:,2) posDest(:,2)];
            Z=[posOrig(:,3) posDest(:,3)];
            
            figure(f(1));
            subplot(length(utypesOrigin),length(utypesDest),cnt);
            
            h2=[];
            for idx=1:length(vals)
                if(m(idx))
                    h2(end+1)=line(X(idx,:),Y(idx,:),Z(idx,:),'Color',colors(idx,:));
                end
            end
            
            %Draw the probe
            link=obj.probe.link;
            s=obj.probe.srcPos;
            d=obj.probe.detPos;
            for iChan = 1:size(link,1)
                iSrc = link.source(iChan);
                iDet = link.detector(iChan);
                
                x = [s(iSrc,1) d(iDet,1)]';
                y = [s(iSrc,2) d(iDet,2)]';
                
                h3 = line(x, y, 'Color', [.7 .7 .7]);
            end
            for i = 1:size(s,1)
                x = s(i,1);
                y = s(i,2);
                text(x, y,['S' num2str(i)], 'FontSize', 14);
            end
            
            for i = 1:size(d,1)
                x = d(i,1);
                y = d(i,2);
                text(x, y,['D' num2str(i)], 'FontSize', 14);
            end
            axis off;
            axis tight;
            title([utypesOrigin{ii} ' --> ' utypesDest{jj}], 'Interpreter','none')
            
            figure(f(2));
            subplot(length(utypesOrigin),length(utypesDest),cnt);
            
            LabelsOrig=strcat(repmat('src-',length(lst),1),num2str(tbl.SourceOrigin(lst)),...
                repmat(':det-',length(lst),1), num2str(tbl.DetectorOrigin(lst))); 
          
            LabelsDet=strcat(repmat('src-',length(lst),1),num2str(tbl.SourceDest(lst)),...
                repmat(':det-',length(lst),1), num2str(tbl.DetectorDest(lst))); 
            [LabelsDet,i]=unique(LabelsDet,'rows');
            [LabelsOrig,j]=unique(LabelsOrig,'rows');
           
          
           vals=reshape(vals(:).*m,length(i),length(j));
           imagesc(vals,[vrange]);
           colorbar;
           set(gca,'XTick',[1:length(i)]);
           set(gca,'YTick',[1:length(j)]);
           
           set(gca,'YTickLabel',{LabelsDet})
           set(gca,'XTickLabel',{LabelsOrig})
            set(gca,'XtickLabelRotation',90);
            title([utypesOrigin{ii} ' --> ' utypesDest{jj}], 'Interpreter','none')
            cnt=cnt+1;
        end
    end
    
    
    