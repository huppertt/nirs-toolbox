function draw( obj, colors, lineStyles, axis_handle )

    % sd pairs
    link = unique( [obj.link.source obj.link.detector], 'rows' );

    s = obj.srcPos;
    d = obj.detPos;
    
    n = size(link, 1);
    
    if nargin < 2 || isempty(colors)
       colors = repmat([0.3 0.5 1], [n 1]);
    elseif size(colors,1) == 1
       colors = repmat(colors, [n 1]);
    end
    
    if nargin < 3 || isempty(lineStyles)
        lineStyles = repmat({'LineStyle', '-', 'LineWidth', 6}, [n 1]);
    elseif size(lineStyles, 1) == 1
        lineStyles = repmat({'LineStyle', '-', 'LineWidth', 6}, [n 1]);
    end
    
    if nargin < 4
        axis_handle = axes();
    end
    
    drawProbe(link, s, d, colors, lineStyles, axis_handle);
    
    labelOptodes( axis_handle, s, d )
    
    rescaleAxes( axis_handle, s, d )
end

function drawProbe(link, s, d, colors, lineStyles, axis_handle)
    axes(axis_handle);
    hold on;
    
    for iChan = 1:size(link,1)
        iSrc = link(iChan,1);
        iDet = link(iChan,2);

        x = [s(iSrc,1) d(iDet,1)]';
        y = [s(iSrc,2) d(iDet,2)]';

        h = line(x, y, 'Color', colors(iChan, :), lineStyles{iChan, :});
        
        set(h,'UserData',[iSrc iDet]);
    end
    
    hold off;
end

function labelOptodes( axis_handle, s, d )
    axes(axis_handle);
    hold on;
     
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
    
    hold off;
end

function rescaleAxes( axis_handle, s, d )
    axes(axis_handle)

    axis equal
    
    p = [s; d];
    
    xmin = min(p(:,1));
    xmax = max(p(:,1));
    
    ymin = min(p(:,2));
    ymax = max(p(:,2));
        
    xl = [xmin xmax];
    yl = [ymin ymax];
    
    xl = 1.2*diff(xl)/2*[-1 1]+mean(xl);
    yl = 1.2*diff(yl)/2*[-1 1]+mean(yl);
    
    axis([xl yl])
    
    axis off
end