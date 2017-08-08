function f = plotSlice(fitobj)
%PLOTSLICE Interactive tool to plot slices of prediction surface
%   plotSlice(fitobj)

%   Copyright 2010-2014 The MathWorks, Inc.


% alpha = 0.05;

% Get information from fit object
[xname,xused,yname,xlims,iscat,catlabels,xfit,xsettings] = getXInfo(fitobj);

plotdata = isscalar(xname);

modname = getString(message('stats:classreg:regr:modelutils:dlg_PredictionSlicePlot'));
slicefig = figure('Units','Normalized','Interruptible','on','Position',[0.05 0.35 0.90 0.5],...
    'NumberTitle','off', 'IntegerHandle','off', ...
    'Name',modname,'Tag','slicefig','ToolBar','none');
fixMenus(slicefig);

ud = fillFigure(slicefig,xlims,iscat,catlabels,fitobj,xsettings,xfit,plotdata,xname,xused,yname);

% Initialize confidence bound menu state
setconf(slicefig,ud.simflag,ud.obsflag,ud.confflag);

set(slicefig,'UserData',ud,'HandleVisibility','callback', ...
    'BusyAction','queue', ...
    'WindowButtonDownFcn',@downFun,...
    'WindowButtonUpFcn',@upFun, ...
    'WindowButtonMotionFcn',@(varargin) motionFun('up',varargin{:}));

% If there are too many x variables, just plot 5 of them to start
if length(xused)>8
    xused(6:end) = [];
    applyAxesSelections(slicefig,xused)
end

if nargout>0
    f = slicefig;
end

% ---------------------------
function ud = fillFigure(slicefig,xlims,iscat,catlabels,fitobj,xsettings,xfit,plotdata,xname,xused,yname,ud)

if nargin<12
    % Define defaults for these settings
    ud.simflag = 1;        % simultaneous confidence intervals
    ud.obsflag = 0;        % not for new observation (for mean)
    ud.confflag = 1;       % show confidence bounds
end

% Create axes and plot fitted lines
slice_axes = makeAxes(slicefig,xused,xlims,iscat,catlabels);
fitline = plotFit(slice_axes,xused,fitobj,xsettings,xfit,iscat,xlims,plotdata,ud);
[ymin,ymax] = updateYAxes(slice_axes);

% Create reference lines (cross-hairs)
reference_line = updateRefLine([],slice_axes,xused,xsettings,ymin,ymax,(ymin+ymax)/2);

% Create UI controls for predictor and response values
[x_field,y_field] = makeUIcontrols(slicefig,xname,xused,yname,iscat,catlabels);

% Fill the rest of users data; save information into figure
ud = makeUserData(xfit,xsettings,xused,fitline,reference_line,slice_axes,...
    iscat,catlabels,fitobj,xlims,x_field,y_field,xname,yname,plotdata,ud);

% Update reference lines
[newy, newyci] = predictionsPlusError(xsettings,fitobj,ud);
updateRefLine(reference_line,slice_axes,xused,xsettings,ymin,ymax,newy);

% Put current predictor and response values into fields
setYField(newy,newyci,ud);
for axnum = 1:length(xused)
    setXField(axnum,xsettings(xused(axnum)),ud);
end

% ----------------------------
function editFun(axnum,tgt,~)
slicefig = ancestor(tgt,'figure');
ud = get(slicefig,'Userdata');
prednum = ud.xused(axnum);

if ~ud.iscat(prednum)
    % Make sure continuous predictor is in range
    val = get(ud.x_field(axnum),'String');
    cx = str2double(val);
    xl = ud.xlims(:,prednum);
    if isnan(cx) || cx<xl(1) || cx>xl(2)
        setXField(axnum,ud.xsettings(prednum),ud);
        warndlg(sprintf('%s',getString(message('stats:classreg:regr:modelutils:dlg_InvalidOrNotInRange',val))), ...
            getString(message('stats:classreg:regr:modelutils:dlg_SlicePlot')),'modal');
        return
    end
end

% Update with new predictor value
ud = updateplot(slicefig, ud, axnum);
set(slicefig,'Userdata',ud);

% ----------------------------
function [y, yci] = predictionsPlusError(x,fitobj,ud)
% Local function for Predicting a response with error bounds.

% Set up data matrix with all variables
in = fitobj.Formula.InModel;
d = fitobj.Variables;
vi = fitobj.VariableInfo;
if isa(d,'dataset')
    d = dataset2table(d);
end
if isa(d,'table')
    d = d(ones(size(x,1),1),:);
    vn = d.Properties.VariableNames;
    cols = find(in);
    for j=1:length(cols)
        vnum = cols(j);
        newdj = d.(vn{vnum});
        if vi.IsCategorical(vnum)
            range = vi.Range{cols(j)};
            for k=1:size(x,1)
                catnum = round(x(k,j));
                if iscell(range)
                    newdj(k) = range(catnum)';
                elseif ischar(range)
                    row = range(catnum,:);
                    newdj(k,:) = ' ';
                    newdj(k,1:length(row)) = row;
                else
                    newdj(k) = range(catnum);
                end
            end
        else
            newdj = x(:,j);
        end
        d.(vn{cols(j)}) = newdj;
    end
else
    d = zeros(size(x,1),length(in));
    d(:,in) = x;
end

% Set up optional arguments
% default: sim false pred curve
args = {};
if ud.simflag
    args(end+(1:2)) = {'Simultaneous' true};
end
if ud.obsflag
    args(end+(1:2)) = {'Prediction' 'observation'};
end

% Predict Response
[y,yci] = predict(fitobj,d,args{:});

% ----------------------------
% Create UI controls on figure to hold predictor and response values
function [x_field,y_field] = makeUIcontrols(slicefig,xname,xused,yname,iscat,catlabels)
fcolor = get(slicefig,'Color');
yfieldp = [.01 .45 .10 .04];
y_field(1) = uicontrol(slicefig,'Style','text','Units','normalized',...
    'Position',yfieldp + [0 .14 0 0],'String','',...
    'ForegroundColor','k','BackgroundColor',fcolor,'Tag','y1');

y_field(2) = uicontrol(slicefig,'Style','text','Units','normalized',...
    'Position',yfieldp + [0 0 0 .04],'String','',...
    'ForegroundColor','k','BackgroundColor',fcolor,'Tag','y2');

uicontrol(slicefig,'Style','Pushbutton','Units','pixels',...
    'Position',[20 10 100 25],'Callback','close','String',getString(message('stats:classreg:regr:modelutils:uicontrol_Close')),'Tag','close');

uicontrol(slicefig,'Style','text','Units','normalized',...
    'Position',yfieldp + [0 0.21 0 .04],'BackgroundColor',fcolor,...
    'ForegroundColor','k','String',yname,'Tag','yname');

n = length(xused);
x_field = zeros(1,n);
for axnum = 1:n
    prednum = xused(axnum);
    xfieldp = [.18 + (axnum-0.5)*.80/n - 0.5*min(.5/n,.15) .09 min(.5/n,.15) .07];
    xtextp  = [.18 + (axnum-0.5)*.80/n - 0.5*min(.5/n,.18) .02 min(.5/n,.18) .05];
    uicontrol(slicefig,'Style','text','Units','normalized',...
        'Position',xtextp,'BackgroundColor',fcolor,...
        'ForegroundColor','k','String',xname{prednum},'Tag',sprintf('xname%d',axnum));
    
    tag = sprintf('xval%d',axnum);
    if iscat(prednum)
        x_field(axnum)  = uicontrol(slicefig,'Style','popup','Units','normalized',...
            'Position',xfieldp,'String',catlabels{prednum},'Tag',tag,...
            'BackgroundColor','white','CallBack',@(varargin) editFun(axnum,varargin{:}));
    else
        x_field(axnum)  = uicontrol(slicefig,'Style','edit','Units','normalized',...
            'Position',xfieldp,'String','','Tag',tag,...
            'BackgroundColor','white','CallBack',@(varargin) editFun(axnum,varargin{:}));
    end
end

% -----------------------
function ud = updateplot(slicefig, ud, axnum)
% Update plot after change in X value or confidence bound setting

slice_axes    = ud.slice_axes;
xused         = ud.xused;
xsettings     = ud.xsettings;
fitline       = ud.fitline;
reference_line= ud.reference_line;
n             = length(slice_axes);

% Get current X values, updating from edit box if necessary
if nargin<3
    axnum = ud.last_axes; % use last axes number if no value specified
end
xrow  = xsettings;
if ~isempty(axnum)
    prednum = xused(axnum);
    cx = getXField(axnum,ud);
    xrow(prednum) = cx;
end
[cy, cyci] = predictionsPlusError(xrow,ud.fitobj,ud);

% If we need to update reference lines, xnum will be non-empty
if ~isempty(axnum)
    xsettings(prednum) = cx;
end

ud.xsettings = xsettings;

for axnum = 1:n
    prednum = xused(axnum);
    xline = getXLine(xsettings,prednum,ud.xfit{prednum},ud.iscat,ud.xlims);
    [yfit, yci] = predictionsPlusError(xline,ud.fitobj,ud);
    
    if ~ud.confflag
        % No conf bounds wanted, so set their y data to NaN so they
        % will not plot but the lines will be around for future use
        yci(:) = NaN;
    end
    set(slice_axes(axnum),'YLimMode','auto');
    set(fitline(1,axnum),'YData',yfit);
    set(fitline(2,axnum),'YData',yci(:,1));
    set(fitline(3,axnum),'YData',yci(:,2));
end

[ymin,ymax] = updateYAxes(slice_axes);
updateRefLine(reference_line,slice_axes,xused,xsettings,ymin,ymax,cy);

ud.last_axes = [];
set(slicefig,'UserData',ud);

setYField(cy,cyci,ud);

% -----------
% helper to locate the axes under the cursor
function axnum = findaxes(fig, allaxes, ~, eventData)

axnum = [];
if matlab.graphics.internal.isGraphicsVersion1()
    h = hittest;
else
    h = eventData.HitObject;
end
if h==fig
    return
end
h = ancestor(h,'axes');
if isempty(h)
    return
end
axnum = find(allaxes==h,1);

% -------------------------
function downFun(varargin)
slicefig = gcbf;
ud = get(slicefig,'Userdata');
set(slicefig,'WindowButtonMotionFcn',@(varargin) motionFun('down',varargin{:}));

axnum = findaxes(slicefig, ud.slice_axes,varargin{:});
if isempty(axnum)
    return
end
ud.last_axes = axnum;
set(slicefig,'Pointer','crosshair');

cp = get(ud.slice_axes(axnum),'CurrentPoint');
cx=cp(1,1);
prednum = ud.xused(axnum);
[xrow,cx] = getXLine(ud.xsettings,prednum,cx,ud.iscat,ud.xlims);
ud.xsettings(prednum) = cx;
[cy, cyci] = predictionsPlusError(xrow,ud.fitobj,ud);

set(slicefig,'Userdata',ud);

setXField(axnum,cx,ud);

set(ud.reference_line(axnum,1),'XData',cx*ones(2,1));
set(ud.reference_line(axnum,2),'YData',[cy cy]);

setYField(cy,cyci,ud);

set(slicefig,'WindowButtonUpFcn',@upFun);

% ----
function motionFun(flag,varargin)
slicefig = gcbf;
ud = get(slicefig,'Userdata');
axnum = findaxes(slicefig, ud.slice_axes,varargin{:});
if isempty(axnum)
    return
end
xrange = get(ud.slice_axes(axnum),'XLim');
newx = getXField(axnum,ud);
maxx = xrange(2);
minx = xrange(1);
n = size(ud.x_field,2);

if isequal(flag,'up') % mouse is up, we are hovering
    % Check for data consistency in each axis
    if n > 1,
        yn = zeros(n,1);
        for idx = 1:n
            y = get(ud.reference_line(idx,2),'Ydata');
            yn(idx) = y(1);
        end
        % If data is inconsistent update all the plots.
        if any(yn~=yn(1))
            upFun(varargin{:});
        end
    end
    % Change cursor to plus sign when mouse is on top of vertical line.
    cursorstate = get(slicefig,'Pointer');
    cp = get(ud.slice_axes(axnum),'CurrentPoint');
    cx = cp(1,1);
    fuzz = 0.02 * (maxx - minx);
    online = cx > newx - fuzz & cx < newx + fuzz;
    if online && strcmp(cursorstate,'arrow'),
        cursorstate = 'crosshair';
    elseif ~online && strcmp(cursorstate,'crosshair'),
        cursorstate = 'arrow';
    end
    set(slicefig,'Pointer',cursorstate);

else % mouse is down, we are dragging
    if ~isequal(ud.last_axes,axnum)
        return;
    end
    cp = get(ud.slice_axes(axnum),'CurrentPoint');
    
    cx=cp(1,1);
    [xrow,cx] = getXLine(ud.xsettings,ud.xused(axnum),cx,ud.iscat,ud.xlims);
    [cy, cyci] = predictionsPlusError(xrow,ud.fitobj,ud);
    
    setXField(axnum,cx,ud);
    
    set(ud.reference_line(axnum,1),'XData', repmat(cx,2,1));
    set(ud.reference_line(axnum,2),'YData',[cy cy]);
    
    setYField(cy,cyci,ud);
end

% ----
function upFun(varargin)
slicefig = gcbf;
set(slicefig,'WindowButtonMotionFcn',@(varargin) motionFun('up',varargin{:}));

ud = get(slicefig,'Userdata');
n = size(ud.x_field,2);
p = get(slicefig,'CurrentPoint');
axnum = floor(1+n*(p(1)-0.18)/.80);

lk = ud.last_axes;
if isempty(lk)
    return
end

xrange = ud.xlims(:,lk);
if axnum < lk
    setXField(lk,xrange(1),ud);
elseif axnum > lk
    setXField(lk,xrange(2),ud);
end

updateplot(slicefig, ud);

% -----------------------------
% Get information about predictors from fit object
function [xname,xused,yname,xlims,iscat,catlabels,xfit,xsettings] = getXInfo(fitobj)
yname = sprintf('%s',getString(message('stats:classreg:regr:modelutils:sprintf_Predicted',fitobj.ResponseName)));
xname = fitobj.PredictorNames;
xused = 1:length(xname);

varinfo = fitobj.VariableInfo;
inmodel = fitobj.Formula.InModel;
catlabels = varinfo.Range(inmodel);
iscat = varinfo.IsCategorical(inmodel);
n = length(iscat);
xlims = ones(2,n);

xfit = cell(1,n);
for prednum=1:n
    range = catlabels{prednum};
    if iscat(prednum)
        if islogical(range)
            labeltext = {'false';'true'};
        elseif isnumeric(range)
            labeltext = num2str(range(:));
        else
            labeltext = char(range);
        end
        catlabels{prednum} = labeltext;
        nlevels = size(labeltext,1);
        xlims(2,prednum) = nlevels;
        xfit{prednum} = 1:nlevels;
    else
        xlims(:,prednum,:) = range(:);
        xfit{prednum} = linspace(range(1),range(2),41);
    end
end

xrange = diff(xlims);
minx = xlims(1,:);

xsettings = minx + xrange/2;
xsettings(iscat) = round(xsettings(iscat));

% ----------------------------
% Create design matrix for predicting over range of predictor prednum using
% current values of other predictors
function [xline,xnew] = getXLine(xsettings,prednum,xfit,iscat,xlims)
xlim = xlims(:,prednum);

xnew = max(xlim(1), min(xlim(2), xfit));
if iscat(prednum)
    xnew = round(xnew);
end

xline = xsettings(ones(length(xnew),1),:);
xline(:,prednum) = xnew;

% ----------------------------
% Get numeric x value from field
function x = getXField(axnum,ud)
if ud.iscat(ud.xused(axnum))
    x = get(ud.x_field(axnum),'Value');
else
    val = get(ud.x_field(axnum),'String');
    x = str2double(val);
end

% ----------------------------
% Set x field from numeric value
function setXField(axnum,x,ud)
if ud.iscat(ud.xused(axnum))
    set(ud.x_field(axnum),'Value',x);
else
    val = num2str(x);
    set(ud.x_field(axnum),'String',val);
end

% ----------------------------
% Set y fields from numeric values
function setYField(y,yci,ud)
set(ud.y_field(1),'String',num2str(double(y)));
set(ud.y_field(2),'String',sprintf('[%g, %g]',yci));

% ----------------------------
% Fix menus on the plotslice menu bar. The toolbar is removed as uicontrols
% are added.
function fixMenus(slicefig)

% Remove menus not appropriate here
menus = findall(slicefig,'type','uimenu');
tags = get(menus,'Tag');
removethese = ismember(tags,{'figMenuTools','figMenuInsert','figMenuView','figMenuEdit'});
delete(menus(removethese));

% Prune file menu
menus = findall(slicefig,'type','uimenu'); % some deleted, find again
menu = findall(menus,'flat','Tag','figMenuFile');
submenus = findall(menus,'flat','Parent',menu);
tags = get(submenus,'Tag');
removethese = ismember(tags, {'figMenuFileExportSetup','figMenuFilePreferences',...
    'figMenuFileSaveWorkspaceAs','figMenuFileImportData','figMenuGenerateCode',...
    'figMenuFileSaveAs','figMenuFileSave','figMenuUpdateFileNew'});
delete(submenus(removethese));
submenus(removethese) = [];
openmenu = findall(submenus,'flat','Tag','figMenuOpen');
delete(openmenu)

% Create menu for controlling confidence bounds
f = uimenu('Label',getString(message('stats:classreg:regr:modelutils:label_Bounds')), 'Position', 2, 'UserData','conf');
uimenu(f,'Label',getString(message('stats:classreg:regr:modelutils:label_Simultaneous')), ...
       'Callback',@doBoundsMenu, 'UserData','simul', 'Tag','boundsSimultaneous');
uimenu(f,'Label',getString(message('stats:classreg:regr:modelutils:label_NonSimultaneous')), ...
       'Callback',@doBoundsMenu, 'UserData','nonsimul', 'Tag','boundsNonsimultaneous');
uimenu(f,'Label',getString(message('stats:classreg:regr:modelutils:label_Curve')), 'Separator','on', ...
       'Callback',@doBoundsMenu, 'UserData','curve', 'Tag','boundsCurve');
uimenu(f,'Label',getString(message('stats:classreg:regr:modelutils:label_Observation')), ...
       'Callback',@doBoundsMenu, 'UserData','observation', 'Tag','boundsObservation');
uimenu(f,'Label',getString(message('stats:classreg:regr:modelutils:label_NoBounds')), 'Separator','on', ...
       'Callback',@doBoundsMenu, 'UserData','none', 'Tag','boundsNone');
   
% Create menu for specifying which predictors are used
f = uimenu('Label',getString(message('stats:classreg:regr:modelutils:label_Predictors')), 'Position', 3, 'UserData','predictors');
uimenu(f,'Label',getString(message('stats:classreg:regr:modelutils:label_Select')),'Callback',@(varargin)SelectPredictors(varargin{:}));

function SelectPredictors(tgt,~)
slicefig = ancestor(tgt,'figure');
ud = get(slicefig,'UserData');

n = length(ud.xsettings);
list = cell(n,1);
for prednum = 1:n
    xname = ud.xname{prednum};
    xval = ud.xsettings(prednum);
    if ud.iscat(prednum)
        labels = ud.catlabels{prednum};
        if iscell(labels)
            xstring = labels{xval};
        elseif isvector(labels)
            xstring = labels(xval);
        else
            xstring = labels(xval,:);
        end
    else
        xstring = num2str(xval);
    end
    list{prednum} = sprintf('%s (%s)',xname,xstring);
end

[selection,ok] = listdlg('ListString',list,'InitialValue',ud.xused,...
                         'PromptString',getString(message('stats:classreg:regr:modelutils:dlg_SlicePlot')));

if ~ok
    return
end

applyAxesSelections(slicefig,selection);

% -----------------------------------------
% Apply selection of a subset of the axes
function applyAxesSelections(slicefig,xused)

clearFigure(slicefig);
ud = get(slicefig,'UserData');

xlims = ud.xlims;
iscat = ud.iscat;
catlabels = ud.catlabels;
fitobj = ud.fitobj;
xsettings = ud.xsettings;
xfit = ud.xfit;

plotdata = ud.plotdata;
xname = ud.xname;
yname = ud.yname;

ud = fillFigure(slicefig,xlims,iscat,catlabels,fitobj,xsettings,xfit,plotdata,xname,xused,yname,ud);
set(slicefig,'UserData',ud);

% ----------------------------
% Clear all figure elements that depend on the data
function clearFigure(slicefig)
ud = get(slicefig,'UserData');
delete(ud.slice_axes)
delete(findall(slicefig,'type','uicontrol'))

% ----------------------------
function doBoundsMenu(tgt,~)
% Update menus to reflect confidence bound settings

menu = tgt;
action = get(menu,'UserData');
slicefig = ancestor(menu,'figure');
ud = get(slicefig,'UserData');

switch(action)
    case 'simul',       ud.simflag = 1;             % simultaneous
    case 'nonsimul',    ud.simflag = 0;             % non-simultaneous
    case 'curve',       ud.obsflag = 0;             % for the mean (fitted line)
    case 'observation', ud.obsflag = 1;             % for a new observation
    case 'none',        ud.confflag = ~ud.confflag; % no confidence intervals
end

% Update menu and bounds
setconf(slicefig,ud.simflag,ud.obsflag,ud.confflag);
updateplot(slicefig, ud);

% ----------------------------
function setconf(slicefig,simul,obs,confflag)
ma = get(findall(slicefig, 'Type','uimenu', 'UserData','conf'), 'Children');
set(ma, 'Checked', 'off');          % uncheck all menu items

% Check item 1 for simultaneous, 2 for non-simultaneous
offon = {'off' 'on'};
set(findobj(ma,'flat', 'Type', 'uimenu', 'UserData', 'simul'),'Checked',offon{1+simul});
set(findobj(ma,'flat', 'Type', 'uimenu', 'UserData', 'nonsimul'),'Checked',offon{2-simul});

% Check item 3 for curve, 4 for new observation
offon = {'off' 'on'};
set(findobj(ma,'flat', 'Type', 'uimenu', 'UserData', 'observation'),'Checked',offon{1+obs});
set(findobj(ma,'flat', 'Type', 'uimenu', 'UserData', 'curve'),'Checked',offon{2-obs});

% Check item 5 if we are omitting confidence bounds.
set(findobj(ma,'flat', 'Type', 'uimenu', 'UserData', 'none'),'Checked',offon{2-confflag});

% ----------------------------
function slice_axes = makeAxes(slicefig,xused,xlims,iscat,catlabels)
n = length(xused);
slice_axes = zeros(n,1);

for axnum = 1:n
    % Create an axis for each input (x) variable
    axisp   = [.18 + (axnum-1)*.80/n .22 .80/n .68];
    
    prednum = xused(axnum);
    slice_axes(axnum) = axes('Parent',slicefig);
    set(slice_axes(axnum),'XLim',xlims(:,prednum)+iscat(prednum)*[-.25;.25],'Box','on','NextPlot','add',...
        'Position',axisp,'GridLineStyle','none');
    if axnum>1
        set(slice_axes(axnum),'Yticklabel',[]);
    end
    if iscat(prednum)
        xlab = catlabels{prednum};
        set(slice_axes(axnum),'XTick',1:length(xlab),'XTickLabel',char(xlab));
    end
end

% ----------------------------
function [ymin,ymax] = updateYAxes(slice_axes)
% Make sure y limits of all axes match; return limits for reference lines

set(slice_axes,'YLimMode','auto'); % force re-calculation

n = length(slice_axes);
yextremes = zeros(n, 2);
for axnum = 1:n
    yextremes(axnum,:) = get(slice_axes(axnum),'YLim');
end

ymin = min(yextremes(:,1));
ymax = max(yextremes(:,2));

set(slice_axes,'YLim',[ymin ymax]);

% ----------------------------
function reference_line = updateRefLine(reference_line,slice_axes,xused,xsettings,ymin,ymax,newy)

if isempty(reference_line)
    % Create Reference Lines
    n = length(slice_axes);
    reference_line = zeros(n,2);
    for axnum = 1:n
        prednum = xused(axnum);
        xlimits = get(slice_axes(axnum),'XLim');
        reference_line(axnum,1) = plot([xsettings(prednum) xsettings(prednum)],[ymin ymax],'--','Parent',slice_axes(axnum));
        reference_line(axnum,2) = plot(xlimits,[newy newy],':','Parent',slice_axes(axnum));
    end
    
    % Reference lines should not cause the axis limits to extend
    set(reference_line(:),'XLimInclude','off','YLimInclude','off');
else
    n = size(reference_line,1);
    for axnum = 1:n
        prednum = xused(axnum);
        set(reference_line(axnum,1),'XData',[xsettings(prednum) xsettings(prednum)],'YData',[ymin ymax]);
        set(reference_line(axnum,2),'YData',[newy newy]); % x limits remain
    end
end

% ----------------------------
function ud = makeUserData(xfit,xsettings,xused,fitline,reference_line,slice_axes,...
                           iscat,catlabels,fitobj,xlims,x_field,y_field,xname,yname,plotdata,ud)
ud.texthandle = [];

ud.xfit = xfit;
ud.xsettings = xsettings;
ud.xused = xused;
ud.fitline = fitline;
ud.reference_line = reference_line;
ud.last_axes = [];
ud.slice_axes = slice_axes;
ud.iscat = iscat;
ud.catlabels = catlabels;
ud.fitobj = fitobj;
ud.x_field = x_field;
ud.y_field = y_field;
ud.xname = xname;
ud.yname = yname;
ud.plotdata = plotdata;

ud.xlims = xlims; % limits of x data, not axis limits

% ----------------------------
% Plotted fitted lines into axes
function fitline = plotFit(slice_axes,xused,fitobj,xsettings,xfit,iscat,xlims,plotdata,ud)

n = length(slice_axes);
fitline = zeros(3,n);

for axnum = 1:n    
    % Show data if requested, which implies just one column in x
    prednum = xused(axnum);
    if plotdata
        ydata = fitobj.Variables.(fitobj.ResponseName);
        if isa(fitobj,'GeneralizedLinearModel') && ...
                strcmpi(fitobj.Distribution.Name,'binomial')
            if size(ydata,2)==2
                n = ydata(:,2);
            else
                n = fitobj.ObservationInfo.BinomSize;
            end
            ydata = ydata(:,1) ./ n;
        end
        line(table2array(fitobj.Variables(:,prednum)),ydata(:,1), 'Linestyle','none', 'Marker','o', 'Parent',slice_axes(axnum));
    end
    
    % Calculate y values for fitted line plot.
    xline = getXLine(xsettings,prednum,xfit{prednum},iscat,xlims);
    [yfit, yci] = predictionsPlusError(xline,fitobj,ud);
    
    % Plot prediction line with confidence intervals
    if iscat(prednum)
        cml = 'go-';
    else
        cml = 'g-';
    end
    fitline(1:3,axnum) = plot(xline(:,prednum),yfit,cml, ...
        xline(:,prednum),yci(:,1),'r:', ...
        xline(:,prednum),yci(:,2),'r:', ...
        'Parent',slice_axes(axnum));
end
