function args = iscoeffwloptimizable(this)
%ISCOEFFWLOPTIMIZABLE True if the object is coeffwloptimizable

%   Copyright 2008 The MathWorks, Inc.

args = [];

% Test if filter is FIR
if ~isfir(this),
    error(message('signal:dfilt:basefilter:iscoeffwloptimizable:FIRRequired'));
end

% Test if designed with IFIR
hmethod = getfmethod(this);
method = hmethod.DesignAlgorithm;
if strcmpi(method,'Interpolated FIR'),
    error(message('signal:dfilt:basefilter:iscoeffwloptimizable:IFIRNotSupported'));
end
    
% Test if the filter was designed by FDESIGN
f = getfdesign(this);
if isempty(f),    
    error(message('signal:dfilt:basefilter:iscoeffwloptimizable:FDESIGNRequired'));
end

% Test supported responses
validresponses = {'Lowpass','Highpass','Halfband','Nyquist'};
if ~any(strcmpi(f.Response,validresponses)),    
    error(message('signal:dfilt:basefilter:iscoeffwloptimizable:InvalidResponse', 'Lowpass', 'Highpass', 'Halfband', 'Nyquist'));    
end

% Verify that the stopband attenuation is defined.
[Fpass, Fstop, Apass, Astop] = minwordlengthspecs(f,this);
if isempty(Astop),
    error(message('signal:dfilt:basefilter:iscoeffwloptimizable:InvalidAstop'));
end

args.fdesignObj = f;
args.Fpass = Fpass;
args.Fstop = Fstop;
args.Apass = Apass;
args.Astop = Astop;

% [EOF]
