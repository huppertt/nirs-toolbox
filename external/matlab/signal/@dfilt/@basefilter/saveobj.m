function s = saveobj(this)
%SAVEOBJ   Save this object.

%   Copyright 1988-2013 The MathWorks, Inc.

s.class   = class(this);
s.version = get(this, 'version');

% Save all of the public properties.
s = setstructfields(s, ...
    savepublicinterface(this));

% Save the reference coefficients.
s = setstructfields(s, ...
    savereferencecoefficients(this));

% Save the metadata.
s = setstructfields(s, ...
    savemetadata(this));

% Save the arithmetic information.
s = setstructfields(s, ...
    savearithmetic(this));

% Save any private data we might need to reproduce the filter.
s = setstructfields(s, ...
    saveprivatedata(this));
  
% Properties added in R2012a and R2012b -----------------------------------
if isprop(this,'FromSysObjFlag') 
  s.FromSysObjFlag = this.FromSysObjFlag;
end

if isprop(this,'FromFilterBuilderFlag') 
  s.FromFilterBuilderFlag = this.FromFilterBuilderFlag;
end

if isprop(this,'ContainedSysObj') && ~isempty(this.ContainedSysObj)
  s.ContainedSysObj = clone(this.ContainedSysObj);
  release(s.ContainedSysObj);
end

if isprop(this,'SupportsNLMethods') 
  s.SupportsNLMethods = this.SupportsNLMethods;
end

% Property added in R2014a ------------------------------------------------
if isprop(this,'FromDesignfilt') 
  s.FromDesignfilt = this.FromDesignfilt;
end


% [EOF]
