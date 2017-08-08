function Hs = tosysobj(this,returnSysObj)
%TOSYSOBJ Convert to a System object

%   Copyright 2012 The MathWorks, Inc.

if ~returnSysObj
  % If returnSysObj is false, then it means that we want to know if the
  % System object conversion is supported for the class at hand. Return a
  % flag as an output instead of throwing an error.
  Hs = false;
  return;
end  

error(message('signal:dfilt:basecatalog:SysobjNotSupported','sysobj',class(this)))
