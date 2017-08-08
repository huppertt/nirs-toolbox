function validatemultiratesysobjstructs(this, method, struct, sysObjFlag)
%DESIGN   

%   Copyright 2011 The MathWorks, Inc.


% If sysObjFlag is true, then check if the requested structure is supported
% by System objects before designing the single rate filter.
if sysObjFlag         
  if any(strcmpi(method,designmethods(this,'iir')))
    error(message('signal:fdesign:basecatalog:IIRStructureNotSupportedBySystemObjects','''SystemObject'''))
  end  
  validStructs = validstructures(this,method,'SystemObject',true);  
  if isempty(intersect(validStructs,struct))
    error(message('signal:fdesign:basecatalog:StructureNotSupportedBySystemObjects',struct,'''SystemObject'''))
  end
end