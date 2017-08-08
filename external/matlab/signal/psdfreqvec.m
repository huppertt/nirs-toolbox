function w = psdfreqvec(varargin)
%This undocumented function may be removed in a future release.
%
%PSDFREQVEC Frequency vector
%   PSDFREQVEC('Npts',NPTS) returns a frequency vector in radians based on
%   the number of points specified in NPTS. The vector returned assumes 2pi
%   periodicity.
%
%   PSDFREQVEC('Fs',FS) specifies the sampling frequency FS in hertz. By
%   default Fs is set to empty indicating normalized frequency.
%   
%   PSDFREQVEC('CenterDC',CENTERDC) specifies a boolean value in CENTERDC
%   which indicates if zero hertz should be in the center of the frequency
%   vector. CENTERDC can be one of these two values [ {false} | true ]. 
%
%   PSDFREQVEC('Range',RANGE) specifies the range of frequency in RANGE.
%   RANGE can be one of the two strings [ {whole} | half ]. Assuming
%   CenterDC=false then:
%       'whole' = [0, 2pi)
%       'half'  = [0, pi] for even NPTS or [0, pi) for odd NPTS
%
%   When CenterDC=true then:
%       'whole' = (-pi, pi] for even NPTs or (-pi, pi) for odd NPTs
%       'half'  = [-pi/2, pi/2] for even* NPTS or (-pi/2, pi/2) for odd NPTS
%
%       *When NPTS is not divisible by 4, then the range is (-pi/2, pi/2).
%
%   When Range='half' the frequency vector has length (NPTS/2+1) if NPTS is
%   even**, and (NPTS+1)/2 if NPTS is odd***.
%
%       **If CenterDc=true and the number of points specified is even is
%       not divisible by 4, then the number of points returned is NPTS/2.
%       This is to avoid frequency points outside the range [-pi/2 pi/2]. 
%
%       ***If CenterDC=true and the number of points NPTS specified is odd
%       and (NPTS+1)/2 is even then the length of the frequency vector is
%       (NPTS-1)/2.

%   Author(s): P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.

% error(nargchk(1,3,nargin,'struct'));

% Define default parameters; use same parameter names as defined in help.
Npts = 1024;
Fs = [];
Range = 'whole';
CenterDC = false;
local_pvparse(varargin{:});

% Compute the frequency grid.
w = frequencygrid(Range,Npts,Fs,CenterDC);

%--------------------------------------------------------------------------
function w = frequencygrid(Range,Npts,Fs,CenterDC)
% Compute the frequency grid.

% Compute the whole frequency range, e.g., [0,2pi) to avoid round off errors.
if isempty(Fs),
    Fs = 2*pi;
end
freq_res = Fs/Npts;
w = freq_res*(0:Npts-1);

% There can still be some minor round off errors in the frequency grid.  
% Fix the known points, i.e., those near pi and 2pi.
Nyq = Fs/2;
half_res = freq_res/2; % half the resolution

% Determine if Npts is odd and calculate half and quarter of Npts.
[isNPTSodd,halfNPTS,ishalfNPTSodd,quarterNPTS] = NPTSinfo(Npts);

if isNPTSodd,
    % Adjust points on either side of Nyquist.
    w(halfNPTS)   = Nyq - half_res;
    w(halfNPTS+1) = Nyq + half_res;
else
    % Make sure we hit Nyquist exactly, i.e., pi or Fs/2 
    w(halfNPTS) = Nyq;
end
w(Npts) = Fs-freq_res;

% Get the right grid based on range, centerdc, etc.
w = finalgrid(w,Npts,Nyq,Range,CenterDC,isNPTSodd,ishalfNPTSodd,halfNPTS,quarterNPTS);

%--------------------------------------------------------------------------
function [isNPTSodd,halfNPTS,ishalfNPTSodd,quarterNPTS] = NPTSinfo(NPTS)
% Determine if we're dealing with even or odd lengths of NPTS, 1/2 NPTS,
% and 1/4 NPTS.

% Determine if Npts is odd.
isNPTSodd = false;
if rem(NPTS,2),
    isNPTSodd = true;
end

% Determine half the number of points.
if isNPTSodd,   halfNPTS = (NPTS+1)/2;  % ODD
else            halfNPTS = (NPTS/2)+1;  % EVEN
end

% Determine if half Npts is odd.
ishalfNPTSodd = false;     
if rem(halfNPTS,2),        
    ishalfNPTSodd = true;  
end

% Determine a quarter of the number of points.
if ishalfNPTSodd,  quarterNPTS = (halfNPTS+1)/2;  % ODD
else               quarterNPTS = (halfNPTS/2)+1;  % EVEN
end

%--------------------------------------------------------------------------
function w = finalgrid(w,Npts,Nyq,Range,CenterDC,isNPTSodd,ishalfNPTSodd,halfNPTS,quarterNPTS)
% Calculate the correct grid based on user specified values for range,
% centerdc, etc.

switch lower(Range)
    case 'whole',
        % Calculated by default.% [0, 2pi)

        if CenterDC,          % (-pi, pi] even or (-pi, pi) odd
            if isNPTSodd,  negEndPt = halfNPTS;
            else           negEndPt = halfNPTS-1;
            end
            w = [-fliplr(w(2:negEndPt)), w(1:halfNPTS)];
        end
        
    case 'half'            
        w = w(1:halfNPTS);      % [0, pi] even or [0, pi) odd
        
        % For even number of points that are not divisible by 4 you get
        % less one point to avoid going outside the [-pi/2 pi/2] range.
        if CenterDC,            % [-pi/2, pi/2] even (-pi/2, pi/2) odd 
            if ishalfNPTSodd,
                negEndPt = quarterNPTS;
            else
                quarterNPTS = quarterNPTS-1; % Avoid going over pi/2
                negEndPt = quarterNPTS;
            end
            w = [-fliplr(w(2:negEndPt)), w(1:quarterNPTS)];
            if ~rem(Npts,4),
                % Make sure we hit pi/2 exactly when Npts is divisible
                % by 4! In this case it's due to roundoff.
                w(end) = Nyq/2;
            end
        end
    otherwise
        error(message('signal:psdfreqvec:InternalError'));
end
w = w(:);  % Return a column vector.

%--------------------------------------------------------------------------
function local_pvparse(varargin)

varnames = evalin('caller','whos');
vnames = {varnames.name};
for m = 1:2:nargin
    indx = find(strncmpi(vnames, varargin{m}, length(varargin{m})));
    switch length(indx)
        case 0
            error(message('signal:psdfreqvec:unknownInput', varargin{ m }));
        case 1
            assignin('caller', vnames{indx}, varargin{m+1});
        otherwise
            error(message('signal:psdfreqvec:ambiguousInput', varargin{ m }));
    end
end

% [EOF]
