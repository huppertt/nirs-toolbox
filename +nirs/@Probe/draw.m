function draw( obj, values, vrange, vcrit, cmap )
    
    % default colormap
    if nargin < 5 || isempty(cmap)
        [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
    end
    
    % no threshold by default
    if nargin < 4
        vcrit = 0;
    end
    
    % default range of values for dispaly
    if nargin == 2
        vmax = ceil( max( abs( values ) ) );
        vrange = [-vmax vmax];
    end
    
    % if we only want to see probe geometry
    if nargin == 1
        cmap    = [0.5 0.5 1];
        vrange  = [0 0];
        values = zeros(size(obj.link.source,1),1);
    end

    % check for scalar vrange or vcrit
    if numel(vrange) == 1
        vrange = [-vrange vrange];
    end
    
    if numel(vcrit) == 1
        vcrit = [-vcrit vcrit];
    end
    
    % mapping of values to colormap index
    vmax = max(abs(vrange));
    z = linspace( -vmax, vmax, size(cmap,1) );

    % threshold colormap
%     lst = z > vcrit(1) & z < vcrit(2);
    
%     [~,i] = min(abs(z));
%     cmap(lst,:) = repmat( cmap(i,:), [sum(lst) 1] );

    lst = values <= vcrit(1) | values >= vcrit(2);

    % loop through channels and draw lines
    link = unique( [obj.link.source obj.link.detector],'rows' );

    s = obj.srcPos;
    d = obj.detPos;

    gca, hold on
    for iChan = 1:size(link,1)

        iSrc = link(iChan,1);
        iDet = link(iChan,2);

        x = [s(iSrc,1) d(iDet,1)]';
        y = [s(iSrc,2) d(iDet,2)]';


        [~, iCol] = min(abs(values(iChan)-z));
        h = line(x,y,'LineStyle','--','LineWidth',4,'Color',cmap(iCol,:));
        
        if lst(iChan) % channel is significant
            set(h,'LineStyle','-','LineWidth',8)
        end
        
        set(h,'UserData',[iSrc iDet])
    end
    
    for i = 1:size(s,1)
        x = s(i,1);
        y = s(i,2);
        
        h = text(x, y,['S' num2str(i)], 'FontSize', 14);
        set(h, 'UserData', ['S' num2str(i)]);
    end
    
    for i = 1:size(d,1)
        x = d(i,1);
        y = d(i,2);
        
        h = text(x, y,['D' num2str(i)], 'FontSize', 14);
        set(h, 'UserData', ['D' num2str(i)]);
    end
    
    axis tight
    axis off

    % colorbar if values were provided
    if nargin > 1
        colorbar
        lst = z > vrange(1) & z < vrange(2);
        colormap(cmap(lst,:))
        caxis(vrange);
    end

    % adjust axes
%     pi = get(gca,'Position');
%     po = get(gca,'OuterPosition');
% 
%     po(3:4) = 1.2*pi(3:4);
%     set(gca,'OuterPosition',po);

    xl = xlim;
    yl = ylim;

    xl = 1.2*diff(xl)/2*[-1 1]+mean(xl);
    yl = 1.2*diff(yl)/2*[-1 1]+mean(yl);

    axis([xl yl])

end