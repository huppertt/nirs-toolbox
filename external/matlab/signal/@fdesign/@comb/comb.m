function this = comb(varargin)
%COMB   Construct a comb filter designer.
%   D = FDESIGN.COMB(COMBTYPE,SPECSTRING,VALUE1,VALUE2,...) constructs a
%   comb filter designer D. Note that D is not the design itself, it only
%   contains the design specifications. In order to design the filter, one
%   needs to invoke the DESIGN method on D. For example (more examples
%   below): 
%   D = fdesign.comb('notch','N,BW',8,.01); 
%   H = design(D); % H is a DFILT
%
%   SPECSTRING is a string that determines what design specifications will
%   be used. There are several possible specifications, a complete list is
%   given below.
%
%   Different specification types may have different design methods
%   available. Use DESIGNMETHODS to get a list of design methods
%   available for a given SPEC: designmethods(D).
%
%   VALUE1, VALUE2, etc. are scalars that provide the value of the
%   corresponding specification. Use get(D,'description') for a description
%   of VALUE1, VALUE2, etc.
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. Moreover, all magnitude specifications are
%   assumed to be in dB.
%
%   COMBTYPE is a string that determines whether the design will correspond
%   to a notching or peaking comb filter. To obtain a notching comb filter
%   design let COMBTYPE = 'notch'. To obtain a peaking comb filter design
%   let COMBTYPE = 'peak'.
%
%   D = FDESIGN.COMB(...,Fs) provides the sampling frequency of the signal
%   to be filtered. Fs must be specified as a scalar trailing all other
%   provided values. For this case, Fs is assumed to be in Hz as are all
%   other frequency values provided. Note that you don't change the
%   specification string in this case. In the example above, if the input
%   signal is sampled at 8 kHz, we can obtain the same filter by specifying
%   the bandwidth in Hz as: D = fdesign.comb('notch','N,BW',8,40,8e3); H =
%   design(D);
%
%   The full list of possible values for SPECSTRING (not case sensitive) is:
%
%      'N,Q' (default)
%      'N,BW'
%      'L,BW,GBW,Nsh'
%
%   where 
%       N     - Filter Order
%       Q     - Peak or Notch Quality Factor  
%       BW    - Peak or Notch Bandwidth
%       GBW   - Gain at which Bandwidth (BW) is measured (dB)
%       L     - Number of Peaks or Notches in the entire [-1 1] normalized
%               frequency interval
%       Nsh   - Shelving Filter Bandwidth 
% 
%   Specification GBW must always be in dB. 
%
%   The comb filter is designed by up-sampling a shelving lowpass (for
%   peaking combs) or highpass (for notching combs) filter design by a
%   factor of L, or N (depending on the SPECSTRING value). The resulting
%   transfer function corresponds to a peak/notch comb filter with L, or N
%   peaks/notches in the Nyquist frequency interval. Designs based on
%   'N,Q', and 'N,BW' specifications use a default unity order shelving
%   filter and a quality factor or bandwidth referenced to a -3 dB
%   gain. Designs based on the 'L,BW,GBW,Nsh' specification allow the
%   choice of an arbitrary value for the shelving filter order Nsh (the
%   default still being unity) as well as an arbitrary gain reference GBW
%   for the peak/notch bandwidth BW (the default value is -3 dB). In this
%   scenario, the number of peaks/notches is chosen by parameter L and the
%   overall comb filter order is equal to Nsh*L. Design parameter Nsh
%   allows the sharpening of the response of the peaks/notches of the comb
%   design.
%
%   D = FDESIGN.COMB uses the  default SPECSTRING ('N,Q') and sets the
%   corresponding values to N, and Q.
%
%   % Example #1 - Design a notching comb filter with 8 notches, and a
%   % notch bandwidth of 0.02 referenced to a -3dB level.
%   d  = fdesign.comb('notch','N,BW',8,0.02)
%   Hd = design(d);
%   fvtool(Hd)
%
%   % Example #2 - Design a notching comb filter with 8 notches, a notch
%   % bandwidth of 0.02 referenced to a -5dB level, and a shelving filter
%   % order of 1.
%   d = fdesign.comb('notch','L,BW,GBW,Nsh',8,0.02,-5,1);
%   Hd = design(d);
%   fvtool(Hd)
%
%   % Example #3 - Design a notching comb filter with 10 notches, a notch
%   % bandwidth of 5 Hz referenced to a -4dB level, a shelving filter
%   % order of 4, and a sampling frequency of 600 Hz.
%   d = fdesign.comb('notch','L,BW,GBW,Nsh',10,5,-4,4,600);
%   Hd = design(d);
%   fvtool(Hd)
%
%   % Example #4 - Design a peaking comb filter with 5 peaks, and a peak
%   % quality factor of 25.
%   d = fdesign.comb('peak','N,Q',5,25);
%   Hd = design(d);
%   fvtool(Hd)
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN, FDESIGN/DESIGNOPTS.

%   Copyright 2008 The MathWorks, Inc.

this = fdesign.comb;

set(this, 'Response', 'Comb Filter');

this.setspecs(varargin{:});

capture(this);

% [EOF]
