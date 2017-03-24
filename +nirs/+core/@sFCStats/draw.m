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
            pval=tbl.pvalue;
            qval=nirs.math.fdr(pval);
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
        case('z')
            values=tbl.Z;
            pval=tbl.pvalue;
            qval=nirs.math.fdr(pval);
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
        otherwise
            warning('type not recognized: using Correlation');
            values=tbl.R;
            pval = tbl.pvalue;
            qval=nirs.math.fdr(pval);
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
            %cmap=flipdim(colormap(autumn(2001)),1);
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
        
        if(~isempty(strfind(s{1},'q')))
            mask = qval < str2double(s{2});
        else
            mask = pval < str2double(s{2});
        end
    end
    I=eye(sqrt(length(mask)));
    mask=mask.*(~I(:));
    
    % range to show
    if nargin < 3 || isempty(vrange)
        vmax    = max(abs(values(:).*mask(:)));
        vrange  = vmax*[-1 1];
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
    
    if(ismember('hyperscan',obj.probe.link.Properties.VariableNames))
        % Draw a hyperscan brain
        p=obj.probe;
        p.link=p.link(ismember(p.link.hyperscan,'A'),:);
        
        S={}; D={};
        for i=1:height(p.link)
            s=['000' num2str(p.link.source(i))];
            S{i}=['Source-' s(end-3:end)];
            d=['000' num2str(p.link.detector(i))];
            D{i}=['Detector-' d(end-3:end)];
        end
        
        p.optodes_registered=p.optodes_registered(ismember(p.optodes.Name,{S{:} D{:}}),:);
        p.optodes=p.optodes(ismember(p.optodes.Name,{S{:} D{:}}),:);
        
        p2=obj.probe;
        p2.link=p2.link(ismember(p2.link.hyperscan,'B'),:);
        
        S={}; D={};
        for i=1:height(p.link)
            s=['000' num2str(p.link.source(i))];
            S{i}=['Source-' s(end-3:end)];
            d=['000' num2str(p.link.detector(i))];
            D{i}=['Detector-' d(end-3:end)];
        end
        
        p2.optodes_registered=p2.optodes_registered(ismember(p2.optodes.Name,{S{:} D{:}}),:);
        p2.optodes=p2.optodes(ismember(p2.optodes.Name,{S{:} D{:}}),:);
        
        for ii=1:length(utypesOrigin)
            for jj=1:length(utypesDest)
                if(strcmp(utypesOrigin(ii),utypesDest(jj)))
                    h1=subplot(2,length(utypesOrigin),ii,'Parent',f(cIdx));
                    s1=p.draw([],[],h1);
                    cb=colorbar(h1,'EastOutside');
                    set(cb,'visible','off');
                    h2=subplot(2,length(utypesOrigin),length(utypesOrigin)+ii,'Parent',f(cIdx));
                    s2=p.draw([],[],h2);
                    cb=colorbar(h2,'EastOutside');
                    set(cb,'visible','off');
                    
                    set(s1,'color',[.5 .5 .5]);
                    set(s2,'color',[.5 .5 .5]);
                    
                    set(h2,'Units','normalized');
                    set(h1,'Units','normalized');
                    
                    if(flip(1)==1);
                        set(h1,'Ydir','reverse');
                        set(h1,'Xdir','normal');
                    else
                        set(h1,'Ydir','normal');
                        set(h1,'Xdir','reverse');
                    end
                    if(flip(2)==1);
                        set(h2,'Ydir','reverse');
                        set(h1,'Xdir','normal');
                    else
                        set(h2,'Ydir','normal');
                        set(h1,'Xdir','reverse');
                    end
                    
                    p.link=p.link(ismember(p.link.type,p.link.type{1}),:);
                    
                    ax=axes(f(cIdx),'Units','normalized','Position',[h1.Position(1) h2.Position(2) h1.Position(3) h1.Position(2)+h1.Position(4)-h2.Position(2)],...
                        'visible','off','Xlim',[-100 100],'Ylim',[-100 100]);
                    hold(ax,'on');
                    axis(ax,'off');
                    hh=getframe(ax);
                    hh=getframe(ax); % For some reason the image is sometimes distorted the first time, but never the 2nd
                    hh=hh.cdata;
                    for i=1:length(s1)
                        s=scatter(h1,s1(i).XData(2),s1(i).YData(2),'filled','r','Sizedata',40);
                        hh2=getframe(ax);
                        [a,b]=find(abs(sum(hh-hh2.cdata,3))>0);
                        i2=p.link.source(i);
                        srcPos(i2,1)=median(b)/size(hh2.cdata,2)*210-105;
                        srcPos(i2,2)=median(a)/size(hh2.cdata,1)*200-100;
                        delete(s);
                        
                        s=scatter(h1,s1(i).XData(1),s1(i).YData(1),'filled','r','Sizedata',40);
                        hh2=getframe(ax);
                        [a,b]=find(abs(sum(hh-hh2.cdata,3))>0);
                        i2=p.link.detector(i);
                        detPos(i2,1)=median(b)/size(hh2.cdata,2)*210-105;
                        detPos(i2,2)=median(a)/size(hh2.cdata,1)*200-100;
                        delete(s);
                        
                    end
                    for i=1:length(s2)
                        s=scatter(h2,s2(i).XData(2),s2(i).YData(2),'filled','r','Sizedata',40);
                        hh2=getframe(ax);
                        [a,b]=find(abs(sum(hh-hh2.cdata,3))>0);
                        i2=p.link.source(i);
                        srcPos2(i2,1)=median(b)/size(hh2.cdata,2)*210-105;
                        srcPos2(i2,2)=median(a)/size(hh2.cdata,1)*200-100;
                        delete(s);
                        
                        s=scatter(h2,s2(i).XData(1),s2(i).YData(1),'filled','r','Sizedata',40);
                        hh2=getframe(ax);
                        [a,b]=find(abs(sum(hh-hh2.cdata,3))>0);
                        i2=p.link.detector(i);
                        detPos2(i2,1)=median(b)/size(hh2.cdata,2)*210-105;
                        detPos2(i2,2)=median(a)/size(hh2.cdata,1)*200-100;
                        delete(s);
                        
                    end
                    srcPos=[srcPos; srcPos2];
                    detPos=[detPos; detPos2];
                    
                    lst=find(strcmp(typesOrigin,utypesOrigin(ii)) & ...
                        strcmp(typesDest,utypesDest(jj)));
                    
                    
                    vals = values(lst);
                    
                    % this mask
                    m = mask(lst);
                    
                    % map to colors
                    idx = bsxfun(@minus, vals', z);
                    [~, idx] = min(abs(idx), [], 1);
                    colors = cmap(idx, :);
                    h2=[];
                    
                    posOrig=(srcPos(tbl.SourceOrigin(lst),:)+...
                        detPos(tbl.DetectorOrigin(lst),:))/2;
                    posDest=(srcPos(tbl.SourceDest(lst),:)+...
                        detPos(tbl.DetectorDest(lst),:))/2;
                    
                    
                    X=[posOrig(:,1) posDest(:,1)];
                    Y=[posOrig(:,2) posDest(:,2)];
                    
                    
                    for idx=1:length(vals)
                        if(m(idx))
                            h2(end+1)=plot(ax,X(idx,:),Y(idx,:),'Color',colors(idx,:));
                        end
                    end
                    set(ax,'YDir','reverse');
                    set(h2,'Linewidth',4)
                    axis(ax,'off');
                    
                    pos=get(ax,'Position');
                    cb=colorbar(ax,'EastOutside');
                    set(ax,'Position',pos);
                    colormap(ax,cmap);
                    caxis(ax,[vrange(1), vrange(2)]);
                    
                    
                    
                end
            end
        end
        
        
        
    else
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
                        posOrig=(obj.probe.srcPos(tbl.SourceOrigin(lst),:)+...
                            obj.probe.detPos(tbl.DetectorOrigin(lst),:))/2;
                        posDest=(obj.probe.srcPos(tbl.SourceDest(lst),:)+...
                            obj.probe.detPos(tbl.DetectorDest(lst),:))/2;
                        
                        
                        X=[posOrig(:,1) posDest(:,1)];
                        Y=[posOrig(:,2) posDest(:,2)];
                        Z=[posOrig(:,3) posDest(:,3)];
                        
                        figure(f(cIdx));
                        subplot(length(utypesOrigin),1,cnt);
                        
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
                        
                    else
                        
                        figure(f(cIdx));
                        subplot(length(utypesOrigin),1,cnt);
                        
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
                    end
                    cnt=cnt+1;
                end
            end
            
        end
    end
    set(f(cIdx),'Name',obj.conditions{cIdx},'NumberTitle','off')
end

