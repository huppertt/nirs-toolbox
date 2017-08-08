function [btype,analog,errStr,msgobj] = iirchk(Wn,varargin)
%IIRCHK  Parameter checking for BUTTER, CHEBY1, CHEBY2, and ELLIP.
%   [btype,analog,errStr] = iirchk(Wn,varargin) returns the 
%   filter type btype (1=lowpass, 2=bandpss, 3=highpass, 4=bandstop)
%   and analog flag analog (0=digital, 1=analog) given the edge
%   frequency Wn (either a one or two element vector) and the
%   optional arguments in varargin.  The variable arguments are 
%   either empty, a one element cell, or a two element cell.
%
%   errStr is empty if no errors are detected; otherwise it contains
%   the error message.  If errStr is not empty, btype and analog
%   are invalid.

%   Copyright 1988-2002 The MathWorks, Inc.

errStr = '';
msgobj = [];

% Define defaults:
analog = 0; % 0=digital, 1=analog
btype = 1;  % 1=lowpass, 2=bandpss, 3=highpass, 4=bandstop

if length(Wn)==1
    btype = 1;  
elseif length(Wn)==2
    btype = 2;
else
    msgobj = message('signal:iirchk:MustBeOneOrTwoElementVector','Wn');
    errStr = getString(msgobj);
    return
end

if length(varargin)>2
    msgobj = message('signal:iirchk:TooManyInputArguments');
    errStr = getString(msgobj);
    return
end

% Interpret and strip off trailing 's' or 'z' argument:
if length(varargin)>0 
    switch lower(varargin{end})
    case 's'
        analog = 1;
        varargin(end) = [];
    case 'z'
        analog = 0;
        varargin(end) = [];
    otherwise
        if length(varargin) > 1
            msgobj = message('signal:iirchk:BadAnalogFlag','z','s');
            errStr = getString(msgobj);
            return
        end
    end
end

% Check for correct Wn limits
if ~analog
   if any(Wn<=0) | any(Wn>=1)
      msgobj = message('signal:iirchk:FreqsMustBeWithinUnitInterval');
      errStr = getString(msgobj);
      return
   end
else
   if any(Wn<=0)
      msgobj = message('signal:iirchk:FreqsMustBePositive');
      errStr = getString(msgobj);
      return
   end
end

% At this point, varargin will either be empty, or contain a single
% band type flag.

if length(varargin)==1   % Interpret filter type argument:
    switch lower(varargin{1})
    case 'low'
        btype = 1;
    case 'bandpass'
        btype = 2;
    case 'high'
        btype = 3;
    case 'stop'
        btype = 4;
    otherwise
        if nargin == 2
            msgobj = message('signal:iirchk:BadOptionString', ...
              'high','stop','low','bandpass','z','s');
            errStr = getString(msgobj);
        else  % nargin == 3
            msgobj = message('signal:iirchk:BadFilterType', ...
              'high','stop','low','bandpass');
            errStr = getString(msgobj);
        end
        return
    end
    switch btype
    case 1
        if length(Wn)~=1
            msgobj = message('signal:iirchk:BadOptionLength','low','Wn',1);
            errStr = getString(msgobj);
            return
        end
    case 2
        if length(Wn)~=2
            msgobj = message('signal:iirchk:BadOptionLength','bandpass','Wn',2);
            errStr = getString(msgobj);
            return
        end
    case 3
        if length(Wn)~=1
            msgobj = message('signal:iirchk:BadOptionLength','high','Wn',1);
            errStr = getString(msgobj);
            return
        end
    case 4
        if length(Wn)~=2
            msgobj = message('signal:iirchk:BadOptionLength','stop','Wn',2);
            errStr = getString(msgobj);
            return
        end
    end
end

