function varargout = design(this, method, varargin)
%DESIGN Design a filter.

%   Copyright 2005-2011 The MathWorks, Inc.

% Get all the design object constructors.
d = getdesignobj(this);

% Make sure that there is a butter field for us to use.
if ~isfield(d, method)
    error(message('signal:fspecs:abstractspec:design:invalidDesign', upper( method )));
end

% Build the object.
d = feval(d.(method));

% Add the System object property if it applies.
supportedStructs = addsysobjdesignopt(d);

% Set any option inputs into the design object.
if nargin > 2
    % First search for a structure of options and set them
    for k = 1:length(varargin),
        if isstruct(varargin{k}),
            set(d,varargin{k});
            varargin = {varargin{1:k-1},varargin{k+1:end}};
            break;
        end
    end
    
    % Set any p-v pairs specified    
    for indx = 1:2:length(varargin)
        if isprop(d, varargin{indx})
            set(d, varargin{indx:indx+1});
        else
            error(message('signal:fspecs:abstractspec:design:invalidOption', varargin{ indx }, upper( method )));
        end
    end
end

% Error out before designing the filter if a System object has been
% requested with a structure that is not supported.
if isprop(d,'SystemObject') && d.SystemObject && ...
    ~any(strcmp(d.FilterStructure,supportedStructs))
    error(message('signal:fspecs:basecatalog:SysObjNotSupported',...
      d.FilterStructure,method,'SystemObject',method,'SystemObject'))
end

% Design the filter.
[varargout{1:nargout}] = design(d, this);

% [EOF]
