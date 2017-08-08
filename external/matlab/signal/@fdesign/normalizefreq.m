function normalizefreq(this)
%NORMALIZEFREQ   Normalize frequency specifications.
%   NORMALIZEFREQ(D) will normalize the frequency specifications in D.
%   The 'NormalizedFrequency' property will be set to true. In addition,
%   all frequency specifications will be normalized by Fs/2 so that they
%   lie between 0 and 1.
%
%   NORMALIZEFREQ(D,BFL) where BFL is either true or false, specifies
%   whether the 'NormalizedFrequency' property will be set to true or
%   false. If not specified, BFL defaults to true. If BFL is set to false,
%   the frequency specifications will be multiplied by Fs/2.
%
%   NORMALIZEFREQ(D,false,Fs) allows for the setting of a new sampling
%   frequency, Fs, when the 'NormalizedFrequency' property is set to false.
%
%   See also FDESIGN, FDESIGN/SETSPECS.

%   Author(s): R. Losada
%   Copyright 2003-2005 The MathWorks, Inc.



% [EOF]
