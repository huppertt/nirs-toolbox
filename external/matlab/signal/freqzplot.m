function freqzplot(h,w,s_in,flag)
%FREQZPLOT Plot frequency response data.
%   FREQZPLOT is obsolete.  FREQZPLOT still works but may be
%   removed in the future. Use FVTOOL instead.
% 
%   See also FREQZ, FVTOOL.

%   Author(s): R. Losada and P. Pacheco 
%   Copyright 1988-2004 The MathWorks, Inc.

warning(message('signal:freqzplot:obsoleteFunction'));

error(nargchk(2,4,nargin,'struct'));

if nargin>3,
    if strcmpi(flag, 'magphase') || strcmpi(flag, 'zerocontphase'),
        phi = h(:,:,2);
        h = h(:,:,1);
    elseif strcmpi(flag, 'zerophase'),
        hr =h;
    end
end        

% Generate defaults
s.xunits  = 'rad/sample'; 
s.yunits  = 'db';
s.plot    = 'both'; % Magnitude and phase
s.fvflag  = 0;
s.yphase  = 'degrees';
if nargin > 2,
    s = parseOpts(s,s_in);
end

% Bring the plot to the foreground
if isfield(s,'ax'),
    ax = s.ax;
    hfig = get(ax,'Parent');
    set(hfig,'CurrentAxes',ax);
else
    
    ax = newplot;
    hfig = get(ax, 'Parent');
    figure(hfig);
    
    hc = get(hfig,'Children');
    for idx = 1:numel(hc)
      newplot(hc(idx))
    end  
    
end

[pd,msg,msgobj] = genplotdata(h,w,s); % Generate the plot data
if ~isempty(msg), error(msgobj); end

switch s.plot,
case 'mag',
    if nargin>3 & strcmpi(flag, 'zerophase'), %#ok
        pd.magh = hr;
        pd.magh = [pd.magh;inf*ones(1,size(pd.magh,2))];
        if length(pd.w)<size(pd.magh,2),
            pd.w = [pd.w;2*pd.w(end)-pd.w(end-1)];
        end
        pd.maglabel = getString(message('signal:freqzplot:Zerophase'));
    end
    plotfresp(ax,pd,'mag');
    
case 'phase',
    if nargin>3,
        pd.phaseh = phi;
        pd.phaseh = [pd.phaseh;inf*ones(1,size(pd.phaseh,2))];
        if length(pd.w)<size(pd.phaseh,1),
            pd.w = [pd.w;2*pd.w(end)-pd.w(end-1)];
        end
        if strcmpi(s.yphase, 'degrees'),
            pd.phaseh = pd.phaseh*180/pi;
        end
        if strcmpi(flag, 'zerocontphase'),
            pd.phaselabel = [getString(message('signal:freqzplot:ContinuousPhase')) ' (' s.yphase ')'];
        end
    end
    plotfresp(ax,pd,'phase');
    
case 'both',
    if nargin>3,
        if size(phi,1) == 1
            phi = phi.';
        end
        pd.phaseh = phi;
        if ~s.fvflag
            pd.phaseh = [pd.phaseh;inf*ones(1,size(pd.phaseh,2))];
        end
        if strcmpi(s.yphase, 'degrees'),
            pd.phaseh = pd.phaseh*180/pi;
        end
    end
    % We plot the phase first to retain the functionality of freqz when hold is on
    ax(2) = subplot(212);
    plotfresp(ax(2),pd,'phase');
    ax(1) = subplot(211);
    plotfresp(ax(1),pd,'mag');
    
    if ishold,
        holdflag = 1;
    else
        holdflag = 0;
    end
    axes(ax(1)); % Bring the plot to the top & make subplot(211) current axis
    
    if ~holdflag,    % Reset the figure so that next plot does not subplot
        set(hfig,'nextplot','replace');       
    end      
end

set(ax,'xgrid','on','ygrid','on','xlim',pd.xlim); 
 
%-----------------------------------------------------------------------------------------
function s = parseOpts(s,s_in)
%PARSEOPTS   Parse optional input params.
%   S is a structure which contains the fields described above plus:
%
%     S.fvflag - flag indicating if a freq. vector was given or nfft was given
%     S.ax     - handle to an axis where the plot will be generated on. (optional)

% Write out all string options
yunits_opts = {'db','linear','squared'};
plot_opts = {'both','mag','phase'};

if ischar(s_in),
	s = charcase(s,s_in,yunits_opts,plot_opts);
	
elseif isstruct(s_in),
	
	s = structcase(s,s_in,yunits_opts,plot_opts);
	
else
    error(message('signal:freqzplot:NeedPlotOpts'));
end

%-------------------------------------------------------------------------------
function s = charcase(s,s_in,yunits_opts,plot_opts)
% This is for backwards compatibility, if a string with freq. units was
% specified as a third input arg. 
indx = find(strncmpi(s_in, yunits_opts, length(s_in)));
if ~isempty(indx),
	s.yunits = yunits_opts{indx};
else
	indx = find(strncmpi(s_in, plot_opts, length(s_in)));
	if ~isempty(indx),
		s.plot = plot_opts{indx};
	else
		% Assume these are user specified x units
		s.xunits = s_in;
	end
end

%-------------------------------------------------------------------------------
function s = structcase(s,s_in,yunits_opts,plot_opts)

if isfield(s_in,'xunits'),
	s.xunits = s_in.xunits;
end

if isfield(s_in,'yunits'),
	s.yunits = s_in.yunits;
end

if isfield(s_in,'plot'),
	s.plot = s_in.plot;
end

if isfield(s_in,'fvflag'),
	s.fvflag = s_in.fvflag;
end

if isfield(s_in,'ax'),
	s.ax = s_in.ax;
end

if isfield(s_in,'yphase'),
    s.yphase = s_in.yphase;
end

% Check for validity of args
if ~ischar(s.xunits),
  error(message('signal:freqzplot:NeedStringFreqUnits'));
end

j = find(strncmpi(s.yunits, yunits_opts, length(s.yunits)));
if isempty(j),
  error(message('signal:freqzplot:InvalidMagnitudeUnits', 'db', 'linear', 'squared'));
end
s.yunits = yunits_opts{j};

k = find(strncmpi(s.plot, plot_opts, length(s.plot)));
if isempty(k),
  error(message('signal:freqzplot:InvalidPlotOpts', 'both', 'mag', 'phase'));
end
s.plot = plot_opts{k};


%-----------------------------------------------------------------------------------------
function plotfresp(ax,pd,type)
switch type
case 'phase'
    data = pd.phaseh;
    lab  = pd.phaselabel;
case 'mag'
    data = pd.magh;
    lab  = pd.maglabel;
    if strcmpi(pd.maglabel, 'zero-phase'),
        lab = 'Amplitude';
    end
end
%axes(ax);
set(ax,'box','on')
line(pd.w,data,'parent',ax);
set(get(ax,'xlabel'),'string',pd.xlab);
set(get(ax,'ylabel'),'string',getTranslatedString('signal:freqzplot',lab));


% [EOF] - freqzplot.m
