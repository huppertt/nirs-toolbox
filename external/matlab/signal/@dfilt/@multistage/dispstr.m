function s = dispstr(this, varargin)
%DISPSTR Coefficient display string for discrete-time filter.
%   DISPSTR(Hd) returns a string representation of the coefficients in the
%   filter. 
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan, J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

if nargin > 1 && isnumeric(varargin{end})
    place = varargin{end};
    varargin(end) = [];
else
    place = [];
end

s = [];

if isempty(place), s = strvcat(get(classhandle(this),'Name'), ' '); end

for indx = 1:nstages(this)
    
    % Build up the "Place" string.
    stagestr = sprintf('%d', indx);
    for jndx = 1:length(place)
        stagestr = sprintf('%d,%s', place(jndx), stagestr);
    end
    stagestr = sprintf('Stage #%s: %s',stagestr,get(classhandle(this.Stage(indx)),'Name'));

    if indx == 1
        spacer = {};
    else
        spacer = {' '};
    end

    % Combine the place string and the info from each section.
    s = strvcat(s, spacer{:}, stagestr, dispstr(this.Stage(indx), varargin{:}, [indx place]));
end

% [EOF]
