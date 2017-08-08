%IIRLP2MB IIR lowpass to M-band frequency transformation.
%   [G, ANUM, ADEN] = IIRLP2MB(H,Wo,Wt,DC) returns transformed filter G for
%   original frequency Wo and target frequency Wt.  The allpass mapping
%   filter is returned in the numerator vector ANUM and the denominator
%   vector ADEN. DC is the string 'pass' or 'stop' to specify which band
%   will be at DC, 'pass' is the default.
%
%   See also FILTERDESIGN/IIRLP2MB.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.



% [EOF]
