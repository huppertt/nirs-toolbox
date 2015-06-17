function draw( obj, values, vcrit, vrange, cmap, axis_handle )
    % draw( )
    % draw( values, vcrit, vrange, cmap, axis_handle )
    
    link = unique( [obj.link.source obj.link.detector], 'rows' );

    s = obj.srcPos;
    d = obj.detPos;
    
    %% no arguments
	if nargin == 1
        axis_handle = axes();
            
        colors  = repmat([0.3 0.5 1], [size(link,1) 1]);
        mask    = logical(ones(size(link,1),1));
        
        drawProbe( axis_handle, link, s, d, colors, mask );

        labelOptodes( axis_handle, s, d );

        rescaleAxes( axis_handle );
    
    %% only axis handle
    elseif nargin == 2 && ishandle(values)
        axis_handle = values; % abusing argument names: values is actually axis_handle
        
        colors  = repmat([0.3 0.5 1], [size(link,1) 1]);
        mask    = logical(ones(size(link,1),1));
        
        drawProbe( axis_handle, link, s, d, colors, mask );

        labelOptodes( axis_handle, s, d );

        rescaleAxes( axis_handle );
    
    %% all arguments
    else
        % no axis handle
        if nargin < 6 || isempty(axis_handle)
            axis_handle = axes();
        end
    
        % masked channels
        if nargin < 3 || isempty(vcrit)
            mask = ones(size(values)) > 0; 
        
        elseif islogical(vcrit)
            mask = vcrit;
            
        elseif isscalar(vcrit)
            mask = values <= -vcrit | values >= vcrit;
            
        else
            mask = values <= vcrit(1) | values >= vcrit(2);
        end
        
        % display range
        if nargin < 4 || isempty(vrange)
            if all( values <= 0 )
                vrange = [min(values) 0];
            elseif all( values >= 0 )
                vrange = [0 max(values)];
            else
                vrange = [-1 1]*max(abs(values));
            end
        end
        
        % colormap
        if nargin < 5 || isempty(cmap)
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
            if all( values <= 0 )
                cmap = cmap(1:ceil(end/2),:);
            elseif all( values >= 0 )
                cmap = cmap(floor(end/2)+1:end,:);
            end
        end
        
        % colors
        z = linspace(vrange(1), vrange(2), size(cmap,1))';
        
        values = min(values, vrange(2));
        values = max(values, vrange(1));
        
        colors = [interp1(z, cmap(:,1), values) ...
            interp1(z, cmap(:,2), values) ...
            interp1(z, cmap(:,3), values)];

        drawProbe( axis_handle, link, s, d, colors, mask );

        labelOptodes( axis_handle, s, d );

        rescaleAxes( axis_handle );

        % colorbar
        c = colorbar;
        a = get(axis_handle, 'Position');

        p = get(c, 'Position');
        p(3) = 0.5*p(3);

        set(c, 'Position', p);

        colormap(cmap);
        caxis(vrange);

        set(axis_handle, 'Position', a);

    end    
end

function drawProbe( axis_handle, link, s, d, colors, mask )
    axes(axis_handle);
    hold on;
    
    for iChan = 1:size(link,1)
        iSrc = link(iChan,1);
        iDet = link(iChan,2);

        x = [s(iSrc,1) d(iDet,1)]';
        y = [s(iSrc,2) d(iDet,2)]';

        h = line(x, y, 'LineStyle', '--', 'LineWidth', 4, 'Color', colors(iChan, :));
        
        if mask(iChan)
            set(h,'LineStyle', '-', 'LineWidth', 8);
        end
        
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

function rescaleAxes( axis_handle )
    axes(axis_handle)

    %axis equal
    %axis tight
    
%     xl = xlim;
%     yl = ylim;
%     
%     xl = 1.2*diff(xl)/2*[-1 1]+mean(xl);
%     yl = 1.2*diff(yl)/2*[-1 1]+mean(yl);
%     
%     axis([xl yl])
    
    axis off
end

