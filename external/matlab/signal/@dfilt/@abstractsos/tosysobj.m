function Hs = tosysobj(this,returnSysObj)
%TOSYSOBJ Convert dfilt SOS structure to System object

%   Copyright 2011 The MathWorks, Inc.

if ~returnSysObj
  % If returnSysObj is false, then it means that we want to know if the
  % System object conversion is supported for the class at hand. Return a
  % flag as an output instead of returning the filter System object.
  Hs = true;
  return;
end  

% Check the number of scale values
checksv(this)

% SOS structures are mapped to a BiquadFilter System object
Hs = dsp.BiquadFilter;

% Call blockparams method to get the structure
s = blockparams(this,'off');

Hs.Structure = s.IIRFiltStruct;

% Make sure the System object coefficients are the refference coefficients
refFilter = reffilter(this);

Hs.SOSMatrixSource = 'Property';
Hs.SOSMatrix = refFilter.sosMatrix;
Hs.ScaleValues = refFilter.ScaleValues;
Hs.OptimizeUnityScaleValues = this.OptimizeScaleValues;

if this.PersistentMemory  
  IC = getinitialconditions(this);
  if isstruct(IC)
    Hs.NumeratorInitialConditions = IC.Num;
    Hs.DenominatorInitialConditions = IC.Den;
  else
    Hs.InitialConditions = IC;
  end
end

if strcmpi(this.Arithmetic, 'fixed')
  Q = get_filterquantizer(this);
  setsysobjfixedpoint(Q,Hs)  
end

% [EOF]