function fdesignhandle = privgetfdesign(this)
%PRIVGETFDESIGN   Get the fdesign handle without copying.
%   PRIVGETFDESIGN should be used when performance is very important.
%   Normally, it is recommended that you use GETFDESIGN as this performs a
%   copy which will protect the FDESIGN metadata from tampering.
%   Unfortunately, this copy is slow when done multiple times.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

fdesignhandle = this.privFDesign;

% [EOF]
