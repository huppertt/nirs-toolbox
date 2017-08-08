function [zh,ph,oh] = zplaneplot(z, p, varargin)
%ZPLANEPLOT Plot Pole/Zero data.
%   ZPLANEPLOT(Z,P) plots the zeros Z and poles P with the unit circle 
%   for reference. ZPLANE(Z,P,AX) puts the plot into the axes specified by
%   the handle AX. 
%
%   ZPLANEPLOT(...,MARKERS) uses the markers defined in the two-element cell 
%   array MARKERS. If Z and/or P are matrices, MARKERS must be specified as a
%   cell array containing elements equal to the number of columns of Z and P,
%   i.e., {'Z1','P1','Z2','P2',...'Zn','Pn'} By default, ZPLANEPLOT represents 
%   a zero with a 'o' and a pole with a 'x'.
%
%   ZPLANEPLOT(...,COLORS) uses the colors defined in the cell array of 
%   three-element RGB vectors COLORS. By default, ZPLANEPLOT uses the 
%   colors specified by the axes ColorOrder property.
%
%   [HZ,HP,Hl] = ZPLANEPLOT(Z,P) returns vectors of handles to the lines and 
%   text objects generated.  HZ is a vector of handles to the zeros lines, 
%   HP is a vector of handles to the poles lines, and Hl is a vector of 
%   handles to the axes / unit circle line and to text objects which are 
%   present when there are multiple zeros or poles.  In case there are no 
%   zeros or no poles, HZ or HP is set to the empty matrix [].
%
%   EXAMPLES:
%    % #1 Custom markers 
%    n = 5; Wn = [100 200]/500;
%    [b,a] = butter(n,Wn);
%    [z,p,k] = tf2zp(b,a);
%    zplaneplot(z,p,{'s','+'});
%    
%    % #2 Custom colors
%    fig = figure; ax = axes;
%    zplaneplot(z,p,{[1 0 0],[.8 .8 .8]});
%
%   % #3 Specify custom markers for two filters (of the same length)
%   n = 20;  a = [1 1 0 0];                 
%   f  = [0 0.4 0.5 1];  f1 = [0 0.65 0.8 1];
%   b  = firpm(n,f,a);    [b1,a1] = eqtflength(b,1);
%   b2  = firpm(n,f1,a);  [b2,a2] = eqtflength(b1,1);
%   [z1,p1,k1] = tf2zp(b1,a1);
%   [z2,p2,k2] = tf2zp(b2,a2);
%   zplaneplot([z1,z2],[p1,p2],{'s','+','d','*'});
%
%  See also ZPLANE, FVTOOL.

%   Copyright 1988-2010 The MathWorks, Inc.

error(nargchk(2,5,nargin,'struct'));

% Parse the optional input arguments.
[ax, markers, color] = parseinputs(varargin{:});

if ~any(imag(z)),
    z = z + 1i*1e-50;
end;
if ~any(imag(p)),
    p = p + 1i*1e-50;
end

% Set some axes properties and return original state.
props = cacheNsetaxesprops(ax);

% Plot the Zeros/Poles using the default axes ColorOrder property for 
% the line colors and the first two elements of the markers cell array.
zh = lclline(z, markers{1}, 7, ax);
ph = lclline(p, markers{2}, 8, ax);

% Use custom markers, if plotting multiple filters.
if length(markers) > 2,
    setProps(zh,ph,'Marker',markers);
end

% Use custom colors, if plotting multiple filters.
if ~isempty(color),
    setProps(zh,ph,'Color',color);
end

% Plot circle.  Calculate the number of points in the polygon based on how
% close the poles and zeroes are to the unit circle.  Do not use less than
% 100 or more than 50000 points.

closest = min(abs(1-[abs(z(:)); abs(p(:))]));

points = 1/closest;
if points < 100
    points = 100;
elseif points > 50000
    points = 50000;
elseif isempty(points)
    points = 100;
end

theta = linspace(0,2*pi,points);
oh = line(cos(theta),sin(theta),'Parent', ax, 'LineStyle', ':');
set(oh, 'HitTest', 'Off'); % Turn hittest off so that the buttondown's will work

set(ax, 'xlimmode', 'auto', 'ylimmode','auto');
get(ax, {'xlim', 'ylim'}); % see g107328

% inline 'axis equal'
units = get(ax,'Units'); set(ax,'Units','Pixels')
apos = get(ax,'Position'); set(ax,'Units',units)
set(ax,'DataAspectRatio',[1 1 1],...
    'PlotBoxAspectRatio',apos([3 4 4]))

set(ax, 'xlimmode', 'auto', 'ylimmode','auto');
get(ax, {'xlim', 'ylim'});

%  zoom out ever so slightly (5%)

if apos(3) < apos(4)
    yl=get(ax,'ylim');
    d=diff(yl);
    yl = [yl(1)-.05*d  yl(2)+.05*d]; 
    set(ax,'ylim',yl);
    xl = get(ax,'xlim');
else
    xl=get(ax,'xlim');
    d=diff(xl);
    xl = [xl(1)-.05*d  xl(2)+.05*d]; 
    set(ax,'xlim',xl); 
    yl = get(ax,'ylim');
end

set(oh,'xdat',[get(oh,'xdat') NaN ...
        xl(1)-diff(xl)*100 xl(2)+diff(xl)*100 NaN 0 0]);
set(oh,'ydat',[get(oh,'ydat') NaN 0 0 NaN ...
        yl(1)-diff(yl)*100 yl(2)+diff(yl)*100]);

oh = [oh drawpznumbers(z, ax)];
oh = [oh drawpznumbers(p, ax)];

% Annotate the axes
set(get(ax,'xlabel'),'string',getString(message('signal:zplaneplot:RealPart')));
set(get(ax,'ylabel'),'string',getString(message('signal:zplaneplot:ImaginaryPart')));

% Reinstate the axes and figure's NextPlot state
resetprops(ax, props);


%-------------------------------------------------------------------
%                       Utility Functions
%-------------------------------------------------------------------
function [ax,markers,color] = parseinputs(varargin)
%
% Outputs:
%   AX      - Handle to the axes which to plot the Pole/Zero plot
%   MARKERS - Cell array of user specified markers for the Poles & Zeros
%   COLOR   - Cell array of user specified RGB color vectors

% Default Markers for Zeros & Poles
markers = {'o','x'};

% Use the default color
color = [];

for n = 1:length(varargin),
    if isequal(length(varargin{n}),1) && ishandle(varargin{n}),
        ax = varargin{n};
        hfig = get(ax, 'Parent');
        visState = get(hfig, 'Visible');
        axes(ax);  %#ok<LAXES>
        set(hfig, 'Visible', visState);
    elseif  iscell(varargin{n}) && ischar(varargin{n}{1}),
        markers = varargin{n}; % Cell array of markers styles 
    elseif  iscell(varargin{n}) && isnumeric(varargin{n}{1})
        color = varargin{n};   % Cell array of color vectors
    end
end

% If no axes handle is specified, create one.
if nargin == 0 || ~exist('ax','var'),
    ax = newplot;
end


%-------------------------------------------------------------------
function props = cacheNsetaxesprops(ax)
%CACHENSETAXESPROPS Set some properties and return original state

% Not sure why original code (PZPLOT local function in ZPLANE)
% had the following code.
% kids = get(ax,'Children');
% for i = 1:length(kids)
%    delete(kids(i));
% end
%set(ax,'xlimmode','auto','ylimmode','auto');


set(ax,'box','on');

% Cache the 'NextPlot' properties
hfig = ancestor(ax,'figure');
props.axnpstate  = get(ax,'NextPlot');
props.fignpstate = get(hfig,'NextPlot');

% equivalent of 'hold on'
set(ax,'NextPlot','Add');
set(hfig,'NextPlot','Add');


%-------------------------------------------------------------------
function resetprops(ax,props)
%RESETPROPS Reinstate the axes and figure's NextPlot state

% Axes specific properties
set(ax,'NextPlot',props.axnpstate);

% Figure specific properties
hfig = ancestor(ax,'figure');
set(hfig,'NextPlot',props.fignpstate);


%-------------------------------------------------------------------
function setProps(zh,ph,propstr,propcell)
% Set Zero/Pole Properties (used for colors and markers)

lenZ = length(zh);
lenP = length(ph);

% Create a vector of inter-laced handles such that we can set
% the appropriate marker style.
h(1:2:2*lenZ) = zh;
h(2:2:2*lenP) = ph;

if isequal(lenZ+lenP,length(propcell)) || lenZ+lenP==0,    
    set(h,{propstr},propcell');       
else
    error(message('signal:zplaneplot:UnMatchedMarkers', propstr));
end


%-------------------------------------------------------------------
function h = lclline(pz, marker, msize, ax)

if isempty(pz)
    pz = NaN;
end
h = line(real(pz), imag(pz), ...
    'LineStyle', 'none', ...
    'Marker', marker, ...
    'MarkerSize', msize, ...
    'Parent', ax);

% [EOF]
