%REORDER   Reorder the sections.
%   REORDER(HD, ORDER) Reorder the sections of the filter using the vector
%   of indices in ORDER.  ORDER does not need to contain all of the indices
%   of the filter, this will remove the missing indices from the filter.  A
%   logical array can also be used to remove sections from the filter, but
%   not to reorder it.  
%
%   HNEW = REORDER(HD, ORDER) If an output is requested a new filter is
%   generated with the new section order and the original filter is not
%   changed.
%
%   REORDER(HD, NUMORDER, DENORDER) Reorders the Numerator and Denominator
%   separately using the vector of indices in NUMORDER and DENORDER
%   respectively.  These vectors must be the same length.
%
%   REORDER(HD, NUMORDER, DENORDER, SVORDER) The ScaleValues can be
%   independently reordered.  If SVORDER is not specified the ScaleValues
%   will be reorder with the Numerator.  The output scale value will always
%   remain on the end when the NUMORDER is used to reorder the scale
%   values.
%
%   REORDER(HD, FILTER_TYPE) where FILTER_TYPE is one of: 'auto','lowpass',
%   'highpass', 'bandpass', or 'bandstop', reorders HD in a way suitable
%   for the given filter type. This mode is intended for fixed-point
%   implementations where the ordering of the sections can have a
%   significant impact in the filter performance.
%
%   Automatic reordering will only take effect when Hd was obtained as a
%   result from a design using FDESIGN. The sections will be automatically
%   reordered depending on the response type of the design (lowpass,
%   highpass, etc.).
%
%   REORDER(HD, DIR_FLAG) If DIR_FLAG is 'UP', the first section will
%   contain the poles closest to the origin, and the last row will contain
%   the poles closest to the unit circle. If DIR_FLAG is 'DOWN', the
%   sections are ordered in the opposite direction. The zeros are always
%   paired with the poles closest to them.
%
%   REORDER(HD, DIR_FLAG, SV) SV is either the string 'poles' or 'zeros'
%   and describes how the ScaleValues should be reordered.  By default the
%   scale values are not reordered when using the DIR_FLAG option.
%
%   % Examples
%   % Use automatic reordering.    
%   Hs = fdesign.lowpass('N,F3db',15, .5);
%   Hd = design(Hs,'butter');
%   reorder(Hd,'auto');
%
%   % Reorder the sections by moving the second section to be in between
%   % the seventh and eighth sections.
%   reorder(Hd, [1 3:7 2 8]);
%   hfvt = fvtool(Hd, 'Analysis', 'coefficients');
%
%   % Remove the third, fourth and seventh sections.
%   Hd1 = reorder(Hd, logical([1 1 0 0 1 1 0 1]));
%   setfilter(hfvt, Hd1);
%
%   % Move the first section to the end and remove the eighth section
%   Hd2 = reorder(Hd, [2:7 1]);
%   setfilter(hfvt, Hd2);
%
%   % Move the numerator and denominator independently.
%   Hd3 = reorder(Hd, [1 3:8 2], [1:8]);
%   setfilter(hfvt, Hd3);
%
%   See also DFILT/CUMSEC, DFILT/SCALE, DFILT/SCALECHECK, FDESIGN.

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

% Help file only

% [EOF]
