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
    vtype='beta';
    
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

if(ismember('hyperscan',obj.probe.link.Properties.VariableNames))
	obj.beta(1:end/2,1:end/2,:) = nan;
    obj.beta(end/2+1:end,end/2+1:end,:) = nan;
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
        case('beta')
            values=tbl.beta;
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
            clabel = 'beta';
        case({'tstat','t'})
            values=tbl.tstat;
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
            clabel = 't-statistic';
        otherwise
            warning('type not recognized: using beta');
            values=tbl.beta;
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
            %cmap=flipdim(colormap(autumn(2001)),1);
    end

    pval=tbl.pvalue;
    qval=nirs.math.BenjaminiHochberg(pval);
    
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
        p = obj.probe;
        
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
                        set(h1,'Xdir','reverse');
                    else
                        set(h1,'Ydir','normal');
                        set(h1,'Xdir','normal');
                    end
                    if(flip(2)==1);
                        set(h2,'Ydir','reverse');
                        set(h2,'Xdir','reverse');
                    else
                        set(h2,'Ydir','normal');
                        set(h2,'Xdir','normal');
                    end
                    
                    p.link=p.link(ismember(p.link.type,p.link.type{1}),:);
                                       
                    lst=find(strcmp(typesOrigin,utypesOrigin(ii)) & ...
                        strcmp(typesDest,utypesDest(jj)));
                                        
                    vals = values(lst);
                    
                    % this mask
                    m = mask(lst);
                    d = tbl(lst,:);
                    
                    ax=axes(f(cIdx),'Units','normalized','Position',[h1.Position(1) h2.Position(2) h1.Position(3) h1.Position(2)+h1.Position(4)-h2.Position(2)],...
                        'visible','off','Xlim',[-100 100],'Ylim',[-100 100]);
                    hold(ax,'on');
                    axis(ax,'off');
                    
                    if any(m)
                    
                        link = p.link;
                        link.X(:,1) = nan;
                        link.Y(:,1) = nan;
                        
                        hh=getframe(ax);
                        hh=getframe(ax); % For some reason the image is sometimes distorted the first time, but never the 2nd
                        hh=hh.cdata;
                        for i=1:length(s1)
                            s=scatter(h1,s1(i).XData(2),s1(i).YData(2),'filled','r','Sizedata',40);
                            hh2=getframe(ax);
                            [a,b]=find(abs(sum(hh-hh2.cdata,3))>0);
                            opt1PosX=median(b)/size(hh2.cdata,2)*210-105;
                            opt1PosY=median(a)/size(hh2.cdata,1)*200-100;
                            delete(s);

                            s=scatter(h1,s1(i).XData(1),s1(i).YData(1),'filled','r','Sizedata',40);
                            hh2=getframe(ax);
                            [a,b]=find(abs(sum(hh-hh2.cdata,3))>0);
                            opt2PosX=median(b)/size(hh2.cdata,2)*210-105;
                            opt2PosY=median(a)/size(hh2.cdata,1)*200-100;
                            delete(s);

                            link.X(i) = (opt1PosX + opt2PosX) / 2;
                            link.Y(i) = (opt1PosY + opt2PosY) / 2;
                        end
                        for i=1:length(s2)
                            s=scatter(h2,s2(i).XData(2),s2(i).YData(2),'filled','r','Sizedata',40);
                            hh2=getframe(ax);
                            [a,b]=find(abs(sum(hh-hh2.cdata,3))>0);
                            opt1PosX=median(b)/size(hh2.cdata,2)*210-105;
                            opt1PosY=median(a)/size(hh2.cdata,1)*200-100;
                            delete(s);

                            s=scatter(h2,s2(i).XData(1),s2(i).YData(1),'filled','r','Sizedata',40);
                            hh2=getframe(ax);
                            [a,b]=find(abs(sum(hh-hh2.cdata,3))>0);
                            opt2PosX=median(b)/size(hh2.cdata,2)*210-105;
                            opt2PosY=median(a)/size(hh2.cdata,1)*200-100;
                            delete(s);
                            
                            link.X(length(s1)+i) = (opt1PosX + opt2PosX) / 2;
                            link.Y(length(s1)+i) = (opt1PosY + opt2PosY) / 2;
                        end

                        % map to colors
                        idx = bsxfun(@minus, vals', z);
                        [~, idx] = min(abs(idx), [], 1);
                        colors = cmap(idx, :);
                        h2=[];

                        % Draw lines with largest magnitude on top
                        [~,sorted_ind] = sort(abs(vals));
                        vals = vals(sorted_ind);
                        m = m(sorted_ind);
                        colors = colors(sorted_ind,:);
                        d = d(sorted_ind,:);
                    end
                                       
                    for idx=1:length(vals)
                        if(m(idx))
                            
                            source_origin = d.SourceOrigin(idx);
                            source_dest = d.SourceDest(idx);
                            detector_origin = d.DetectorOrigin(idx);
                            detector_dest = d.DetectorDest(idx);
                            
                            origin_ind = find(link.source == source_origin & link.detector == detector_origin);
                            dest_ind = find(link.source == source_dest & link.detector == detector_dest);
                            
                            X_origin = link.X( origin_ind );
                            Y_origin = link.Y( origin_ind );
                            X_dest = link.X( dest_ind );
                            Y_dest = link.Y( dest_ind );

                            h2(end+1)=plot(ax,[X_origin X_dest],[Y_origin Y_dest],'Color',colors(idx,:));
                            
                        end
                    end

                    set(ax,'YDir','reverse');
                    set(h2,'Linewidth',4)
                    axis(ax,'off');
                    
                    pos=get(ax,'Position');
                    cb=colorbar(ax,'EastOutside');
                    set(ax,'Position',pos);
                    if ~any(m)
                        colormap(ax,[0 0 0]);
                        set(cb,'ytick',[0 1],'yticklabel',{'','n.s.'})
                        ylabel(cb,' ');
                    else
                        colormap(ax,cmap);
                        ylabel(cb,clabel);
                    end
                    caxis(ax,[vrange(1), vrange(2)]);
                    
                    if strcmp( utypesOrigin{ii} , utypesDest{ii} )
                        title( utypesOrigin(ii) );
                    else
                        title( [utypesOrigin{ii} ' - ' utypesDest{ii} ] );
                    end
                    
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
    supertitle(f(cIdx),obj.conditions{cIdx});
end

