function Hd = cheby2(this, varargin)
%CHEBY2 Chebyshev Type II digital filter design.

%   Copyright 1988-2011 The MathWorks, Inc.

N = this.FilterOrder;
Fs = this.Fs;
Fc = this.F3dB;
Ast = this.Astop;

% Compute analog frequency
if this.NormalizedFrequency,
    Wc = tan(pi*Fc/2);
else
    Wc = tan(pi*Fc/Fs);
end

% Find corresponding analog stopband-edge frequency
Wst = Wc*cosh(1/N*acosh(sqrt(10^(Ast/10)-1)));

% Convert analog stopband-edge frequency to digital
Fst = 2*atan(Wst)/pi;

% Construct a new fdesign object with the converted specs
Hdes2 = fdesign.lowpass('N,Fst,Ast',N,Fst,Ast);

% If System object has been requested, do not design a System object here.
% Let the privdesigngateway method convert the dfilt object to a System
% object.
[varargin sysObjFlag] = parsesysobj(Hdes2,'design',varargin{:});
  
Hd = cheby2(Hdes2,varargin{:});

% Set the SystemObject design option property to the value in sysObjFlag
d = getfmethod(Hd);
if isprop(d,'SystemObject')
  d.SystemObject = sysObjFlag;
end

% Error out if a System object has been requested with a structure that is
% not supported.
supportedStructs = getvalidsysobjstructures(d);
method = getdesignmethod(Hd);
if isprop(d,'SystemObject') && d.SystemObject && ...
    ~any(strcmp(d.FilterStructure,supportedStructs))
    error(message('signal:fspecs:basecatalog:SysObjNotSupported',...
      d.FilterStructure,method,'SystemObject',method,'SystemObject'))
end

% [EOF]
