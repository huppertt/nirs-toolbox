function dopts = designoptions(this, method, varargin)
%DESIGNOPTIONS Return the design options.

%   Copyright 2004-2011 The MathWorks, Inc.


if nargin < 2,
    error(message('signal:fdesign:abstracttypewspecs:designoptions:notEnoughInputs', inputname( 1 )));
end

idx = strcmpi(varargin,'signalonly');
sigOnlyFlag = ~isempty(idx);
varargin(idx) = [];

% Parse the SystemObject input
sysObjFlag = validatedesignoptionssysobjinput(this,varargin{:});

% Only try to call DESIGNOPTIONS on the specs if the method is valid.
if isdesignmethod(this, method)
    if sigOnlyFlag
      dopts = designoptions(this.CurrentSpecs, method, true); 
      
      if isfield(dopts,'DensityFactor')
        dopts = rmfield(dopts,{'DensityFactor','DefaultDensityFactor'});
      end
      
      if isfield(dopts,'SystemObject')
        dopts = rmfield(dopts,{'SystemObject','DefaultSystemObject'});
      end
            
    else
      dopts = designoptions(this.CurrentSpecs, method);   
    end
    % If SystemObject parameter was passed as an input with a value of true,
    % then remove structures that are not supported by System objects.
    if sysObjFlag
      hf = feval(getdesignobj(getcurrentspecs(this),method));
      s = getvalidsysobjstructures(hf);      
      if isempty(s)
        % Error out if none of the structures is supported by System
        % objects
        error(message('signal:fdesign:basecatalog:NoSupportedSysObj',method,'''SystemObject'''))
      end
      dopts.FilterStructure = s;
    end    
else
    error(message('signal:fdesign:abstracttypewspecs:designoptions:invalidMethod', upper( method ), this.Specification));
end

% [EOF]
