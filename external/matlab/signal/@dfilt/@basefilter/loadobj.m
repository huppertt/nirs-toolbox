function this = loadobj(s)
%LOADOBJ   Load this object.

%   Copyright 1988-2013 The MathWorks, Inc.

this = feval(s.class);

% Fix old versions.
if isstruct(s) && ~isfield(s, 'version')
    s.version.number      = 0;
    s.version.description = 'R14';
end

loadreferencecoefficients(this, s);

% Set the arithmetic before the public interface in case subclasses have
% properties in the public interface which are actually in the arithmetic
% (the filter quantizers).
loadarithmetic(this, s);

loadpublicinterface(this, s);

% Load the private data after all public properties to make sure setting
% the public property doesn't overwrite the private settings.
loadprivatedata(this, s);

% Load the metadata last because it affects nothing.
loadmetadata(this, s);

% Load System object related properties added in R2012a, R2012b
if isstruct(s)
  if isfield(s,'FromSysObjFlag')
    this.FromSysObjFlag = s.FromSysObjFlag;
  end
  if isfield(s,'FromFilterBuilderFlag')
    this.FromFilterBuilderFlag = s.FromFilterBuilderFlag;
  end  
  if isfield(s,'ContainedSysObj') && ~isempty(s.ContainedSysObj)
    this.ContainedSysObj = clone(s.ContainedSysObj);
    release(this.ContainedSysObj)
  end
  if isfield(s,'SupportsNLMethods')
    this.SupportsNLMethods = s.SupportsNLMethods;
  end   
  % designfilt related property added in R2014a
  if isfield(s,'FromDesignfilt')
    this.FromDesignfilt = s.FromDesignfilt;
  end  
else
  if isprop(this,'FromSysObjFlag')
    this.FromSysObjFlag = s.FromSysObjFlag;
  end
  if isprop(s,'FromFilterBuilderFlag')
    this.FromFilterBuilderFlag = s.FromFilterBuilderFlag;
  end    
  if isprop(this,'ContainedSysObj') && ~isempty(s.ContainedSysObj)
    this.ContainedSysObj = clone(s.ContainedSysObj);
    release(this.ContainedSysObj)
  end
  if isprop(this,'SupportsNLMethods')
    this.SupportsNLMethods = s.SupportsNLMethods;
  end  
  % designfilt related property added in R2014a
  if isprop(this,'FromDesignfilt')
    this.FromDesignfilt = s.FromDesignfilt;
  end  
  
end
  
end
% [EOF]
