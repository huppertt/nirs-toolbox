function varargout=draw( obj, colors, lineStyles, axis_handle,norescale)
    %% draw - Plots the probe geometry.
    % 
    % Args:
    %     colors      - (optional) n x 3 array of colors [R, G, B] for each channel
    %     lineStyles  - (optional) 2D cell array with each row containing
    %                   arguments for the 'line' functions (e.g. {'LineWidth',6})
    %     axis_handle - (optional) handle to axis to the plot to
        
    % sd pairs


    if(nargin<5)
        norescale=false;
    end
    [~,unique_pairs_idx]=unique(obj.link(:,[1,2]),'rows');
    link=obj.link(sort(unique_pairs_idx),1:2);


    s = obj.srcPos;
    d = obj.detPos;
    
    n = height(link);
    
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
    
    h=drawProbe(link, s, d, colors, lineStyles, axis_handle);
    
    labelOptodes( axis_handle, s, d );
    
    if(~norescale)
        rescaleAxes( axis_handle, s, d );
    end
    if(nargout==1)
        varargout{1}=h;
    end
    
end

function h=drawProbe(link, s, d, colors, lineStyles, axis_handle)
%    axes(axis_handle);
    hold(axis_handle,'on');
    h = [];
    ss=scatter(s(:,1),s(:,2),'filled','MarkerFaceColor','r');
    sd=scatter(d(:,1),d(:,2),'filled','MarkerFaceColor','b');
    set(ss,'Tag','SourceOptodes');
    set(sd,'Tag','DetectorOptodes');
            

    for iChan = 1:height(link)
        if iscell(link.source(iChan))
            iSrc = link.source{iChan};
            iDet = link.detector{iChan};
        else
            iSrc = link.source(iChan);
            iDet = link.detector(iChan);
        end

        for j = 1:length(iSrc)
            x = [s(iSrc(j),1) d(iDet(j),1)]';
            y = [s(iSrc(j),2) d(iDet(j),2)]';

            h(iChan) = line(x, y,'parent',axis_handle, 'Color', colors(iChan, :), lineStyles{iChan, :});
        end
        
        set(h(iChan),'UserData',[iSrc iDet]);
    end
    
    hold(axis_handle,'off');
end

function labelOptodes( axis_handle, s, d )
 %   axes(axis_handle);
    hold(axis_handle,'on');
     
    for i = 1:size(s,1)
        x = s(i,1);
        y = s(i,2);
        
        h = text(x, y,['S' num2str(i)],'parent',axis_handle, 'FontSize', 14);
        set(h, 'UserData', ['S' num2str(i)]);
    end
    
    for i = 1:size(d,1)
        x = d(i,1);
        y = d(i,2);
        
        h = text(x, y,['D' num2str(i)],'parent',axis_handle, 'FontSize', 14);
        set(h, 'UserData', ['D' num2str(i)]);
    end
    
    hold(axis_handle,'off');
end

function rescaleAxes( axis_handle, s, d )
   % axes(axis_handle)

    axis(axis_handle,'equal');
    
    p = [s; d];
    
    xmin = min(p(:,1));
    xmax = max(p(:,1));
    
    ymin = min(p(:,2));
    ymax = max(p(:,2));
        
    xl = [xmin xmax];
    yl = [ymin ymax];
    
    if ~all( [diff(xl) diff(yl)] > 0 )
       xl = xlim;
       yl = ylim;
    end
    
    xl = 1.2*diff(xl)/2*[-1 1]+mean(xl);
    yl = 1.2*diff(yl)/2*[-1 1]+mean(yl);
    
    axis(axis_handle,[xl yl])
    
    axis(axis_handle,'off')
end