function varargout = fxptdesign(this, method, varargin)
%FXPTDESIGN   Design a minimum word length filter.
%   FXPTDESIGN(D, M, VARARGIN) Design the filter using the method in the string
%   M on the specs D.  

%   Copyright 2008 The MathWorks, Inc.

% Test if Fixed-Point Designer is installed
if ~isfixptinstalled,
     error(message('signal:fdesign:lowpass:fxptdesign:fixptTbxRq'));
end

% Error out if not FIR
if isempty(find(designmethods(this,'fir'),method)),
     error(message('signal:fdesign:lowpass:fxptdesign:FIRonly'));
end

% Minimum wordlength design (i.e. Smart default for CoeffWordLength)
Hd = design(this,method,varargin{:});
Hd = minwordfir(this,Hd,'noiseShaping',ns,newargs{:});

% [EOF]
