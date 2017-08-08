%VALIDSTRUCTURES   Return the valid structures
%   V = VALIDSTRUCTURES(D) returns a structure, V, that contains a list of
%   all the valid filter structures for each design method available for
%   the object D.
%
%   V = VALIDSTRUCTURES(D, METHOD) returns a cell array, V, containing the
%   valid filter structures for the object D and the design method, METHOD.
%
%   V = VALIDSTRUCTURES(D, ..., 'SystemObject', FLAG) limits the list of
%   filter structures to those that are supported by System objects when
%   FLAG is set to true. Lists all the filter structures supported by
%   DFILT/MFILT objects when FLAG is set to false. Omitting the
%   'SystemObject', FLAG input pair is equivalent to setting the FLAG to
%   false (DSP System Toolbox required).
%
%   % Example
%   h = fdesign.lowpass;
%   validstructures(h, 'equiripple')
%
%   See also FDESIGN/DESIGNMETHODS and FDESIGN/DESIGNOPTS.

% Copyright 2004-2011 The MathWorks, Inc.
