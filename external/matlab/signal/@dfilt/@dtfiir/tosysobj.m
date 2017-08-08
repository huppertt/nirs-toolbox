function Hs = tosysobj(this,returnSysObj)
%TOSYSOBJ Convert dfilt IIR structure to System object

%   Copyright 2012 The MathWorks, Inc.

if ~returnSysObj
  % If returnSysObj is false, then it means that we want to know if the
  % System object conversion is supported for the class at hand. Return a
  % flag as an output instead of returning the filter System object.
  Hs = true;
  return;
end  

% IIR structures are mapped to an IIRFilter System object
Hs = dsp.IIRFilter;
Hs.Structure = getFilterStructure(this);

% Make sure the System object coefficients are the refference coefficients
refFilter = reffilter(this);

Hs.Numerator = refFilter.Numerator;
Hs.Denominator = refFilter.Denominator;

if this.PersistentMemory  
  IC = getinitialconditions(this);  
  if ~isempty(IC)
    if isPropertyActive(Hs,'InitialConditions')
      Hs.InitialConditions = IC;
    else
      Hs.NumeratorInitialConditions = IC.Num;
      Hs.DenominatorInitialConditions = IC.Den;
    end
  end
end

if strcmpi(this.Arithmetic, 'fixed')
  Q = get_filterquantizer(this);
  setsysobjfixedpoint(Q,Hs)  
end

function s = getFilterStructure(this)
s = this.FilterStructure;
s = strrep(s,'-',' ');
s = strrep(s,'F','f');

% [EOF]