function [pd,msg,msgobj] = genplotdata(h,w,s)
%GENPLOTDATA Generate frequency response plotting data.
%
%   [PD,MSG] = GENPLOTDATA(H,W,S) generates the plotting data for 
%   the frequency response H computed at the frequencies specified in 
%   the vector W (in rad/sample).  S specifies additional plotting 
%   information that can be altered for different plotting options.  
%
%   PD is a structure which contains some of the following fields:
%   PD.W          - Frequency data
%   PD.XLAB       - Frequency axis (x-axis) label
%   PD.XLIM       - Frequency axis (x-axis) limits
%   PD.MAGH       - Magnitude data (possibly in dB)
%   PD.MAGLABEL   - Magnitude label
%   PD.PHASEH     - Unwrapped phase data (possibly in radians)
%   PD.PHASELABEL - Phase label

%   Author(s): R. Losada and P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.

% If h is a vector, make it a column
if ndims(h) == 2 && min(size(h)) == 1
  h = h(:);
end

% Initialize outputs
 pd = [];
 msg = '';
 msgobj = [];
 
% Force w to be a column vector.
if ~any(size(w)==1),
   msgobj = message('signal:genplotdata:FreqInputMustBeVector','W');
   msg = getString(msgobj);
   return;   
end
pd.w = w(:);

changed_freq = 0; % This flag is used when s.plot is set to 'both'. If the freq vector,
                  % pd.w, is changed, the flag will be set to 1, this way we won't
                  % change pd.w again when generating the phase data.
                      
% Generate the appropriate data according to the desired plot type

% Cell array of the standard magnitude strings (used for the Ylabels)
magunitstrs = getmagunitstrs;

if any(strncmp(s.plot,{'mag','both'}, length(s.plot))),
    % Generate the magnitude data
    pd.magh     = abs(h);
    pd.maglabel = magunitstrs{1};
    
    % Convert to Magnitude (dB)
    if strncmpi(s.yunits, 'db', length(s.yunits)),
        pd.magh     = convert2db(pd.magh);
        pd.maglabel =  magunitstrs{2};
        
    % Convert to Magnitude Squared    
    elseif strncmpi(s.yunits, 'squared', length(s.yunits)),
        pd.magh     =  convert2sq(pd.magh);
        pd.maglabel =  magunitstrs{3};
    end
    
    % Make sure you show the nyquist or 2*nyquist since freqz doesn't return this value
    if ~s.fvflag,
        pd.w = add_freq_point(pd.w);
        changed_freq = 1;
        pd.magh = [pd.magh;inf*ones(1,size(pd.magh,2))]; 
    end
end

if any(strncmp(s.plot, {'phase','both'}, length(s.plot))),
    % Generate the phase data
    pd.phaseh = unwrap(angle(h));
    if strcmpi(s.yphase, 'degrees'),
        pd.phaseh = pd.phaseh*180/pi;
    end
    pd.phaselabel = ['Phase (' s.yphase ')'];
    
    % Make sure you show the nyquist or 2*nyquist since freqz doesn't return this value
    if ~s.fvflag,
        if ~changed_freq,
            pd.w = add_freq_point(pd.w);    
        end
        pd.phaseh = [pd.phaseh;inf*ones(1,size(pd.phaseh,2))]; 
    end
end


% Generate the correct frequency units and label     

% Cell array of the standard frequency units strings (used for the Xlabels)
frequnitstrs = getfrequnitstrs;

switch lower(s.xunits),
case 'rad/sample',
    pd.xlab = frequnitstrs{1};
    pd.w    = pd.w./pi; % Scale by pi in the plot
case 'hz',
    pd.xlab = frequnitstrs{2}; 
case 'khz',
    pd.xlab = frequnitstrs{3};
case 'mhz',
    pd.xlab = frequnitstrs{4}; 
case 'ghz',
    pd.xlab = frequnitstrs{5};
otherwise
    pd.xlab = s.xunits;
end

% Set x-axis range to be exactly the start and end points of w
pd.xlim = [pd.w(1) pd.w(end)]; 


%-----------------------------------------------------------------------------------------
function w_out = add_freq_point(w)
%ADD_FREQ_POINT   adds an extra frequency point to the frequency vetor.
%   To be used when the frequency vector flag of FREQZ is not set to
%   make sure the entire Nyquist interval is shown in the plot even if
%   the Nyquist point is not included.
%
%   We made this a local function even though it is a one-liner because
%   it is called in two different places.

w_out = [w;2*w(end)-w(end-1)];

% [EOF]
