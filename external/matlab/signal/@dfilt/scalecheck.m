function scalecheck(Hd)
%SCALECHECK   Check the scaling of an SOS filter.
%   S = SCALECHECK(Hd,PNORM) for DF1SOS and DF2TSOS filters returns a row
%   vector S with that computes the p-norm of the filter from its input to
%   the output of each second-order section. Note that this computation
%   does not include the trailing scale value Hd.ScaleValue(end).
%
%   PNORM can be either frequency-domain norms: 'L1', 'L2', 'Linf' or
%   discrete-time-domain norms: 'l1', 'l2', 'linf'. Note that the L2-norm
%   of a filter is equal to its l2-norm (Parseval's theorem), but this is
%   not true for other norms.
%
%   For DF2SOS and DF1TSOS, S is a row vector containing the p-norm from
%   the filter input to the input of the recursive part of each
%   second-order section. This corresponds to the input to the multipliers
%   in these structures and are the points where overflow should be
%   avoided. If Hd has non-trivial scale values, i.e. if not all scale
%   values are equal to one, S is a two-row matrix containing in the
%   second row the p-norm from the input of the filter to the input of each
%   scale value between the sections. Note that for these structures, the
%   last numerator and the trailing scale value are not included when
%   checking the scale.
%
%   For a given p-norm, an optimally scaled filter will have partial norms
%   equal to one, so S will contain all ones.
%
%     Example: Check the Linf-norm scaling of a filter
%     Hs = fdesign.lowpass; % Create a filter design specifications object.
%     Hd = design(Hs,'ellip');       % Design an elliptic SOS filter
%     scale(Hd,'Linf');
%     scalecheck(Hd,'Linf')
%
%   See also DFILT/SCALE, DFILT/NORM, and DFILT/REORDER

%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.



% [EOF]
