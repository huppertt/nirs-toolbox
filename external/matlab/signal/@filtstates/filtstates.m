function varargout = filtstates(varargin)
%FILTSTATES   Filter States object.
%   H = FILTSTATES.DFIIR(NUMSTATES,DENSTATES) constructs an object and sets
%   its 'Numerator' and 'Denominator' properties to NUMSTATES and DENSTATES
%   respectively.  
%
%   Notice that the DSP System Toolbox, along with the Fixed-Point
%   Designer, enables single precision floating-point and fixed-point
%   support for the Numerator and Denominator states.
%
%   The following methods are available for the DFIIR object (type
%   help filtstates/METHOD to get help on a specific method - e.g. help
%   filtstates/double):
%
%   filtstates/double - Convert a FILTSTATES object to a double vector.
%   filtstates/single - Convert a FILTSTATES object to a single vector.
%
%   For more information, enter doc filtstates at the MATLAB command line.
    
%   Author(s): P. Costa
%   Copyright 1988-2012 The MathWorks, Inc.

% Help for instantiating a FILTSTATES object.

error(message('signal:filtstates:filtstates:Package'));

% [EOF]














