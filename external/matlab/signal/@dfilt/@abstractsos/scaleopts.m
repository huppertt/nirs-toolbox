function opts = scaleopts(this)
%SCALEOPTS   Create an options object for second-order section scaling.
%   OPTS = SCALEOPTS(Hd) create an options object OPTS that contains scaling
%   options for second-order section scaling. Arithmetic specific values are set
%   according to the settings in Hd. The OPTS object can be passed in as an
%   argument to the SCALE method.
%
%   See also DFILT/SCALE, DFILT/REORDER, DFILT/SCALECHECK, DFILT/CUMSEC,
%   DFILT/NORM.

%   Copyright 2003-2009 The MathWorks, Inc.

% Construct default opts object
opts = fdopts.sosscaling;

% Set arithmetic specific defaults
scaleopts(this.filterquantizer,this,opts);

% [EOF]
