function varargout = thisdesign(this, method, varargin)
%DESIGN   

%   Copyright 2005-2011 The MathWorks, Inc.

% Parse the inputs for the filterstructure for fir designs
struct = [];
if ~any(strcmpi(method,designmethods(this,'iir'))),
    [struct, varargin] = parsestruct(this, 'firdecim', method, varargin{:});
end

% If SystemObject has been passed as an input, remove it and cache its
% value. We do not want to convert the single rate filter design to a
% System object. We will convert the final mfilt object. Validate that the
% requested structure is supported by System objects.
[varargin sysObjFlag] = parsesysobj(this,'design',varargin{:});
validatemultiratesysobjstructs(this, method, struct, sysObjFlag)

n = nargout;
if n == 0
    n = 1;
end

% Call the DESIGN method of the contained FDesign object.
if strcmpi(method, 'multisection')
    [varargout{1:n}] = design(this.CurrentFDesign, 'multisection', this.DecimationFactor);
else
    [varargout{1:n}] = design(this.CurrentFDesign, method, varargin{:});
end

Hm = varargout{1};

% If the filter is not already in an interpolation object, use one.
if ~isa(Hm, 'mfilt.abstractmultirate')
    
    M = get(this, 'DecimationFactor');
    
    % Get the coefficients from the filter.
    b = tf(Hm);
    
    % Cache the FMETHOD object handle.
    hfmethod = Hm.getfmethod;
    
    % Set the SystemObject property, if it exists, to the cached sysObjFlag
    % value
    if isprop(hfmethod,'SystemObject')
      hfmethod.SystemObject = sysObjFlag;
    end

    Hm = feval(['mfilt.' struct], M, b);
    Hm.setfmethod(hfmethod);
            
else
  % If it is an mfilt object, then get the fmethod object and add the
  % System object design option property
  fm = getfmethod(Hm);
  if ~isprop(fm,'SystemObject')
    p = schema.prop(fm, 'SystemObject', 'bool');
    set(p, ...
    'FactoryValue',          false, ...
    'AccessFlags.AbortSet',  'Off', ...
    'AccessFlags.Serialize', 'Off', ...
    'AccessFlags.Copy',      'Off');   
  end
  fm.SystemObject = sysObjFlag;
end

varargout{1} = Hm;
% [EOF]
