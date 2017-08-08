classdef BoxPlotListenerManager < matlab.graphics.internal.SaveGraphicsVersion2Only 
    
    properties
        FigureHandle
        AxesHandle
        BoxParentHandle
        LabelAxis
    end
    
    properties(SetAccess = private, GetAccess = private, Transient)
        Listeners
    end
    
    properties(SetAccess = private, GetAccess = private)
        % must be the last property defined
        % it triggers the creation of the listeners during deserialization
        ListenersAdded = false;
    end
    
    methods
        function this = BoxPlotListenerManager(f, ax, boxparent, labelAxis)
            this = this@matlab.graphics.internal.SaveGraphicsVersion2Only;
            bplm = getappdata(boxparent,'BoxPlotListenerManager');
            if isa(bplm,'BoxPlotListenerManager')
                this = bplm;
            end
            this.FigureHandle = f;
            this.AxesHandle = ax;
            this.BoxParentHandle = boxparent;
            this.LabelAxis = labelAxis;
            this.ListenersAdded = true; % Setting this property creates the listeners.
            setappdata(boxparent,'BoxPlotListenerManager',this);
        end
        
        function delete(this)
            internal.stats.BoxPlotListenerManager.removeListeners(this.Listeners);
        end

        function set.ListenersAdded(this, val)
            if ~this.ListenersAdded && val
                this.ListenersAdded = true;
                this.Listeners = internal.stats.BoxPlotListenerManager.makeListeners(this.FigureHandle, this.AxesHandle, this.LabelAxis);
            end
        end
    end
    
    methods (Static)
        
        function destroy(boxparent) 
            bplm = getappdata(boxparent,'BoxPlotListenerManager');
            if isa(bplm,'BoxPlotListenerManager')
                setappdata(boxparent,'BoxPlotListenerManager',[]);
            end
        end
            
        %
        %----------------------------
        % Set factor-axis label positions and visibility, and axes position.
        % Shift labels to their correct spot, and adjust the axes size to
        %  make sure the labels are visible. Permit the labels to be painted, now
        %  that they are in the right spots. Leave invisible those labels not
        %  required by the labelverbosity setting, or those that are off the edge
        %  of the axis.
        function repositionLabels(f,ax)

            boxparent = getappdata(ax,'boxplothandle');

            % Guaranteed to be called only when a boxplot hggroup is in the axes, but
            % doublecheck anyway.
            if isempty(boxparent) || ~ishghandle(boxparent)
                return; % Not a boxplot.
            end
            plotType = getappdata(boxparent,'plottype');
            if ~strcmp(plotType,'boxplot')
                return; % Not a boxplot.
            end

            % Fetch the appdata.
            labelPtsPos = getappdata(boxparent,'labelptspos');
            labelDatLoc = getappdata(boxparent,'labeldatloc');
            labelHandles = getappdata(boxparent,'labelhandles');
            columnPtsPosition = getappdata(boxparent,'columnptsposition');
            labelAxis = getappdata(boxparent,'labelaxis');
            displayLabel = getappdata(boxparent,'displaylabel');

            switch(labelAxis)
              case 'x'
                dat = 1;
                edge = 2;
                datspan = 3;
                span = 4;
                labellimtag = 'XLim';
                labelscaletag = 'XScale';
                labeldirtag = 'XDir';
                axislabel = 'XLabel';
                axisticklabel = 'XTickLabel';
                axislimoffset = 1;
              case 'y'
                edge = 1;
                dat  = 2;
                span = 3;
                datspan = 4;
                labellimtag='YLim';
                labelscaletag = 'YScale';
                labeldirtag = 'YDir';
                axislabel = 'YLabel';
                axisticklabel = 'YTickLabel';
                axislimoffset = 3;
              otherwise
                % Something isn't set up right.
                return; % Abort from this axes with no warning.
            end

            % Check that each handle is still valid.
            validHandles = ishghandle(labelHandles);
            if ~any(validHandles)
                return; % No work to do, try another axes.
            end

            % Skip invalid handles, in case certain labels were deleted.
            if ~all(validHandles)
                labelHandles = labelHandles(validHandles);
                labelPtsPos = labelPtsPos(validHandles);
                labelDatLoc = labelDatLoc(validHandles);
                displayLabel = displayLabel(validHandles);
            end

            % Get axis label handle.
            haxlabel = get(ax,axislabel);
            % Get tick label text.
            ticklabel = get(ax,axisticklabel);

            % Get the current axis position values in axis units
            axUnits = get(ax,'Units');
            dftlUnits = get(f,'DefaultAxesUnits');
            li = get(ax,'LooseInset');
            lid = get(f,'DefaultAxesLooseInset');
            op = get(ax,'OuterPosition');
            if ~isequal(axUnits,dftlUnits)
                lid = hgconvertunits(f,lid,dftlUnits,axUnits,get(ax,'Parent'));
            end

            % Convert the axes position values to points.
            liPts = hgconvertunits(f,li, axUnits,'points',get(ax,'Parent'));
            lidPts = hgconvertunits(f,lid, axUnits,'points',get(ax,'Parent'));
            opPts = hgconvertunits(f,op, axUnits,'points',get(ax,'Parent'));

            % Determine how much room is available for tick labels, assuming the axes
            % were squished to zero width.
            maxAllowedWidthPts = opPts(span)-lidPts(edge)-lidPts(span);

            % Determine how many factor labels we can fit.
            if columnPtsPosition(1)<columnPtsPosition(end)
                % Inline labels.
                columnPosShifted = -columnPtsPosition(1)+[columnPtsPosition(2:end),0];
            else
                % Stacked labels, or just one factor.
                columnPosShifted = -columnPtsPosition;
            end
            okcols = columnPosShifted<=maxAllowedWidthPts;
            if ~isempty(okcols) && dat==1  % always show at least one label on x axis
                okcols(1) = 1;
            end
            numDisplayedFactors = find(okcols,1,'last');

            % Determine how much room the labels require, and which specific labels
            % will be displayed.  Labels are omitted if there is not enough room on the
            % axes, beginning with the minor labels.
            if isempty(numDisplayedFactors) || ~isempty(ticklabel)
                % All labels clipped.
                labelClipped = true(size(labelPtsPos));
                if isempty(ticklabel)
                    actualWidthPts = 0; % if no factors displayed or empty tick labels
                else
                    % Need to get measurements from current tick labels instead of
                    % using ones from text measured originally
                    if isequal(labelAxis,'y')
                        sampletext = ticklabel; % to get max width of all labels
                        ind = 3;
                    else
                        sampletext = '42'; % any single-line text will do
                        ind = 4;
                    end
                    h = text(0,0,sampletext,'Parent',ax,'Visible','off', ...
                             'Interpreter','none','Units','points');
                    ext = get(h,'Extent');
                    actualWidthPts = ext(ind);
                    delete(h);
                end
            elseif numDisplayedFactors==length(columnPosShifted)
                % No labels clipped.
                labelClipped = false(size(labelPtsPos));
                actualWidthPts = columnPosShifted(end);
            else
                % Some, but not all, factor tick labels don't fit and will be clipped.
                if columnPtsPosition(1)<columnPtsPosition(end)
                    % Inline labels.
                    labelClipped = labelPtsPos>columnPtsPosition(numDisplayedFactors);
                    % Slide the labels closer to the axis, to fill the gap left by
                    % omitted labels.
                    labelPtsPos = labelPtsPos...
                        +columnPosShifted(end)...
                        -columnPosShifted(numDisplayedFactors);
                else
                    % Stacked labels.
                    labelClipped = labelPtsPos<columnPtsPosition(numDisplayedFactors);
                end
                actualWidthPts = columnPosShifted(numDisplayedFactors);
            end

            % Adjust axes to accommodate the factor tick labels and any axis label.
            % Set the other three edges to their default.
            newLiPts(dat) = lidPts(dat);
            newLiPts([span datspan]) = liPts([span datspan]);
            newLiPts(edge) = lidPts(edge)+actualWidthPts;

            % Convert back to original axes units.
            newLi = hgconvertunits(f,newLiPts, ...
                                   'points',get(ax,'Units'),get(ax,'Parent'));

            % If the tick labels take up too much room, squish axis to 0 width but do
            % not let its width go negative.  Do this in the original units, to avoid
            % late roundoff.
            if newLi(edge)+newLi(span)>op(span)
                newLi(edge) = max(0, op(span)-newLi(span));
            end

            % Set the margin needed around axes; the axes Position property
            % will adjust itself.
            set(ax,'LooseInset',newLi);

            % Get the updated axes position.
            p = get(ax,'Position');
            pPts = hgconvertunits(f,p, ...
                                  get(ax,'Units'),'points',get(ax,'Parent'));

            % For non-warped axes (as in "axis square"), recalculate another way
            if isequal(get(ax,'WarpToFill'),'off')
                xl = get(ax,'XLim');
                yl = get(ax,'YLim');
                
                % Use text to get coordinate (in points) of southwest corner
                t1 = text(xl(1),yl(1),'42','Visible','off');
                set(t1,'Units','points');
                pSW = get(t1,'Position');
                delete(t1);
                
                % Same for northeast corner
                t1 = text(xl(2),yl(2),'42','Visible','off');
                set(t1,'Units','points');
                pNE = get(t1,'Position');
                delete(t1);
                
                % Re-create position; we only care about the last two elements
                % Use min/max/abs in case one or more directions are reversed
                pPts = [min(pSW(1),pNE(1)),...
                        max(pSW(2),pNE(2)),...
                        abs(pNE(1)-pSW(1)),...
                        abs(pNE(2)-pSW(2))];
            end

            % Convert label tick locations from data to points units.

            % Get the factor axis settings.
            scale = get(ax,labelscaletag);
            dir = get(ax,labeldirtag);
            lim = get(ax,labellimtag);
            if strcmp(scale,'log') && any(lim==0)
                % Find the actual limits by getting the hidden deprecated RenderLimits.
                oldstate = warning('off',...
                                   'MATLAB:HandleGraphics:NonfunctionalProperty:RenderLimits');
                renderlimits = get(ax,'RenderLimits');
                warning(oldstate);
                lim = renderlimits([axislimoffset,axislimoffset+1]);
            end

            % Map the label tick locations from data units into normalized units.
            switch scale
              case 'linear'
                labelDatLocNorm = (labelDatLoc-lim(1))/(lim(2)-lim(1));
              case 'log'
                labelDatLocNorm = ...
                    (log(labelDatLoc)-log(lim(1)))/(log(lim(2))-log(lim(1)));
              otherwise
                error(message('stats:boxplot:BadScale'));
            end
            % Flip the direction, if requested.
            switch dir
              case 'normal' %do nothing
              case 'reverse'
                labelDatLocNorm = 1-labelDatLocNorm;
              otherwise
                error(message('stats:boxplot:BadDir'));
            end
            % Find which labels are outside the axis limits, so they will be made
            % invisible.
            labelOutOfRange = labelDatLocNorm<0 | labelDatLocNorm>1 ...
                | isnan(labelDatLocNorm);

            % Map the normalized units into points units.
            axisLengthPts = pPts(datspan);
            labelDatLocPts = labelDatLocNorm*axisLengthPts;

            % Fill a cell array with the label positions.
            labelpos = cell(size(labelPtsPos));
            onepos = zeros(1,3);
            onepos(3) = -.1;  % Slight negative Z, so datatip appears in front.
            for j=1:length(labelpos)
                onepos(dat) = labelDatLocPts(j);
                onepos(edge) = labelPtsPos(j);
                labelpos{j} = onepos;
            end

            % Set the label positions.  Be sure to specify points units, as
            % the units may have changed (for example, during printing).
            set(labelHandles,'Units','points',{'Position'},labelpos);

            % Make labels within the limits visible, and those outside invisible.
            % displayLabel makes some labels always invisible, based on the
            % labelverbosity setting.
            vis = displayLabel & ~labelOutOfRange & ~labelClipped;
            visoptions = {'off';'on'};
            visvalue = visoptions(vis+1);
            set(labelHandles,{'Visible'},visvalue);

            % Position the axis label
            % Use slight negative Z, so datatip appears in front.
            axlabelPtsposition(3) = -.1;
            axlabelPtsposition(edge) = -(actualWidthPts+1);
            axlabelPtsposition(dat) = pPts(datspan)/2;
            set(haxlabel,'Position',axlabelPtsposition);

        end

        function listeners = makeListeners(f, ax, labelAxis) 
            % Fire callback if window is resized.
            list1 = addlistener(f,'SizeChanged',@(src,evt) internal.stats.BoxPlotListenerManager.repositionLabels(f,ax));

            % Fire callback if certain axes properties change, eg axis limits due to
            % zooming.
            props = [{'DataAspectRatioMode' 'DataAspectRatio' 'WarpToFill' ...
                      'PlotBoxAspectRatio' 'XLim' 'YLim' 'Position'}, ...
                     strcat(upper(labelAxis),{'Dir' 'Scale' 'TickLabel'})];
            list2 = addlistener(ax, props, 'PostSet',@(src,evt) internal.stats.BoxPlotListenerManager.repositionLabels(f,ax));
            
            % If units change, temporarily ignore other changes      
            list3 = addlistener(ax, 'Units', 'PreSet',@(src,evt) enableListenerCallback(list2,false));
            list4 = addlistener(ax, 'Units', 'PostSet',@(src,evt) enableListenerCallback(list2,true));
             
            list5 = addlistener(ax, 'ObjectBeingDestroyed', @(src,evt)delete(list1));
            listeners = {list1 list2 list3 list4 list5};
        end

        function removeListeners(listeners)
            for j=1:numel(listeners)
                listj = listeners{j};
                if isa(listj,'handle.listener') || isa(listj,'property.listener')
                    delete(listj)
                end
            end
        end
    end
end

function enableListenerCallback(listener,isEnabled)
if matlab.graphics.internal.isGraphicsVersion1
    if isEnabled
        onoffFlag  = 'on';
    else
        onoffFlag = 'off';
    end
    set(listener,'Enabled',onoffFlag)
else
    listener.Enabled = isEnabled;
end
end
