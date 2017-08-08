function Hs = tosysobj(this,returnSysObj)
%TOSYSOBJ Convert dfilt Lattice AR structure to System object

%   Copyright 2012 The MathWorks, Inc.

if ~returnSysObj
  % If returnSysObj is false, then it means that we want to know if the
  % System object conversion is supported for the class at hand. Return a
  % flag as an output instead of returning the filter System object.
  Hs = true;
  return;
end  

% Lattice AR structures are mapped to an Allpole System object
Hs = dsp.AllpoleFilter('Structure','Lattice AR');

% Make sure the System object coefficients are the refference coefficients
refFilter = reffilter(this);

Hs.ReflectionCoefficients = refFilter.Lattice;

if this.PersistentMemory  
  IC = getinitialconditions(this);  
  if ~isempty(IC)
    Hs.InitialConditions = IC;
  end
end

if strcmpi(this.Arithmetic, 'fixed')
  Q = get_filterquantizer(this);
  setsysobjfixedpoint(Q,Hs)  
end

% [EOF]