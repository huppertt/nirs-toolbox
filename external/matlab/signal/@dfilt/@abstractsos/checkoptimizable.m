function checkoptimizable(this, hTar)
%   OUT = CHECKOPTIMIZABLE(ARGS) Check if optimizations can be carried on

%   Copyright 2009 The MathWorks, Inc.

if any(this.ScaleValues==0) && strcmpi(hTar.OptimizeZeros,'on'),
   warning(message('signal:dfilt:abstractsos:checkoptimizable:NullScaleValues'));
    hTar.OptimizeZeros = 'off';
end


% [EOF]
