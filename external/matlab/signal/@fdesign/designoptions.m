%DESIGNOPTIONS  Show all design options available for a particular design
%   S = DESIGNOPTIONS(D, METHOD) returns a structure, S, with all design
%   options available for an object D using a particular design method,
%   METHOD.
%   
%   S = DESIGNOPTIONS(D, METHOD, 'SystemObject', FLAG) limits the list of
%   filter structures to those that are supported by System objects when
%   FLAG is set to true. Lists all the filter structures supported by
%   DFILT/MFILT objects when FLAG is set to false. Omitting the
%   'SystemObject', FLAG input pair is equivalent to setting the FLAG to
%   false (DSP System Toolbox required).
%
%   % Example
%   h = fdesign.lowpass;
%   designoptions(h, 'equiripple')
%
%   See also FDESIGN/DESIGNMETHODS and FDESIGN/VALIDSTRUCTURES.

%   Copyright 2006-2011 The MathWorks, Inc.



% [EOF]
