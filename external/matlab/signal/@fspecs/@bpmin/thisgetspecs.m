function specs = thisgetspecs(this)
%THISGETSPECS   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fstop1 = this.Fstop1;
specs.Fpass1 = this.Fpass1;
specs.Fpass2 = this.Fpass2;
specs.Fstop2 = this.Fstop2;
specs.Astop1 = this.Astop1;
specs.Apass  = this.Apass;
specs.Astop2 = this.Astop2;

% specs.fpass = [this.Fpass1 this.Fpass2];
% specs.fstop = [0 this.Fstop1 this.Fstop2 1];
% specs.apass = this.Apass;
% specs.astop = [this.Astop1 this.Astop2];

% specs.stopband.frequency = [0 this.Fstop1];
% specs.stopband.astop     = this.Astop1;
% 
% specs.passband.frequency = [this.Fpass1 this.Fpass2];
% specs.passband.apass     = this.Apass;
% 
% specs.stopband(2).frequency = [this.Fstop2 1];
% specs.stopband(2).astop     = this.Astop2;

% [EOF]
