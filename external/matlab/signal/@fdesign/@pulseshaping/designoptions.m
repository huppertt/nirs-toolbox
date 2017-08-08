function dopts = designoptions(this, method, varargin)
%DESIGNOPTIONS Return the design options.

%   Copyright 2008-2011 The MathWorks, Inc.

% Parse the SystemObject input
sysObjFlag = validatedesignoptionssysobjinput(this,varargin{:});

dopts = designoptions(this.PulseShapeObj, method);

% If SystemObject parameter was passed as an input with a value of true,
% then remove structures that are not supported by System objects.
if sysObjFlag
  hf = getfmethod(this,method);
  if ismultistage(hf)
    s = [];
  else
    supportedStructs = getsysobjsupportedstructs(hf);  
    [s k ~] = intersect(dopts.FilterStructure,supportedStructs);        
  end
  if isempty(s)
    % Error out if none of the structures is supported by System objects
    error(message('signal:fdesign:basecatalog:NoSupportedSysObj',method,'''SystemObject'''))
  end  
  dopts.FilterStructure = dopts.FilterStructure(sort(k)); 
end

% [EOF]
