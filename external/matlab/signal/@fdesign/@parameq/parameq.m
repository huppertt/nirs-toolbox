function this = parameq(varargin)
%PARAMEQ   Construct a parametric equalizer filter designer.
%   D = FDESIGN.PARAMEQ(SPECSTRING,VALUE1,VALUE2,...) constructs a
%   parametric equalizer filter designer D. Note that D is not the design
%   itself, it only contains the design specifications. In order to design
%   the filter, one needs to invoke the DESIGN method on D.
%   For example (more examples below):
%   D = fdesign.parameq('N,F0,BW,Gref,G0,GBW',4,0.3,0.2,0,10,7);
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
%   corresponding specification. In the example above, this means that N =
%   4, F0 = 0.3, BW = 0.2, Gref = 0, G0 = 10, and GBW = 7. Use 
%   get(D,'description') for a description of VALUE1, VALUE2, etc.
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. Moreover, all magnitude specifications are
%   assumed to be in dB.
%
%   D = FDESIGN.PARAMEQ(...,Fs) provides the sampling frequency of the
%   signal to be filtered. Fs must be specified as a scalar trailing the
%   other numerical values provided. For this case, Fs is assumed to be in
%   Hz as are all other frequency values provided. Note that you don't
%   change the specification string in this case. In the example above, if
%   the input signal is sampled at 8 kHz, we can obtain the same filter by
%   specifying the frequencies in Hz as:
%   D = fdesign.parameq('N,F0,BW,Gref,G0,GBW',4,1200,800,0,10,7,8e3);
%   H = design(D);
%
%   Unlike most filter design objects (fdesign), specifications for
%   parametric equalizers must be specified in dB. Also, the gains are the
%   actual levels (in dB) and can be either positive or negative.
%   Parametric equalizers can either boost (G0 > Gref) or cut (Gref > G0)
%   the input signal. In addition to that, the gains must satisfy either   
%   G0 > Gp > GBW > Gst > Gref (boost) or
%   G0 < Gp < GBW < Gst < Gref (cut) 
%   for any of such values that is specified. Note that you do not always
%   specify all gains. The gains specified depend on SPECSTRING. Also, the
%   bandwidths must satisfy:
%   BWst > BW > BWp
%
%   The full list of possible values for SPECSTRING (not case sensitive) is:
%
%      'F0,BW,BWp,Gref,G0,GBW,Gp' (min. order default)
%      'F0,BW,BWst,Gref,G0,GBW,Gst'
%      'F0,BW,BWp,Gref,G0,GBW,Gp,Gst'
%      'N,F0,BW,Gref,G0,GBW'
%      'N,F0,BW,Gref,G0,GBW,Gp'
%      'N,F0,BW,Gref,G0,GBW,Gst'
%      'N,F0,BW,Gref,G0,GBW,Gp,Gst'
%      'N,F0,Fc,Qa,G0'
%      'N,F0,Fc,S,G0'
%      'N,F0,Qa,Gref,G0'
%      'N,Flow,Fhigh,Gref,G0,GBW'
%      'N,Flow,Fhigh,Gref,G0,GBW,Gp'
%      'N,Flow,Fhigh,Gref,G0,GBW,Gst'
%      'N,Flow,Fhigh,Gref,G0,GBW,Gp,Gst' 
%
%   where 
%       F0    - Center Frequency
%       Fc    - Cutoff Frequency (for shelving lowpass/highpass filters)
%       BW    - Bandwidth
%       Qa    - Quality Factor (for audio-related specifications)
%       S     - Slope parameter (for shelving lowpass/highpass filters)
%       BWp   - Passband Bandwidth
%       BWst  - Stopband Bandwidth
%       Gref  - Reference Gain (dB)
%       G0    - Center Frequency Gain (dB)
%       GBW   - Gain at which Bandwidth (BW) is measured (dB)
%       Gp    - Passband Gain (dB)
%       Gst   - Stopband Gain (dB)
%       N     - Filter Order
%       Flow  - Lower Frequency at Gain GBW
%       Fhigh - Higher Frequency at Gain GBW
%
%   D = FDESIGN.PARAMEQ(F0,BW,BWp,Gref,G0,GBW,Gp) uses the  default
%   SPECSTRING ('F0,BW,BWp,Gref,G0,GBW,Gp') and sets the
%   corresponding values to F0, BW, etc.
%
%   % Example #1 - Design a parametric equalizer filter that boosts by 6 dB
%   d  = fdesign.parameq('F0,BW,BWp,Gref,G0,GBW,Gp',...
%      0.4,0.3,0.2,0,6,6+10*log10(0.5),5);
%   Hd = design(d);
%   fvtool(Hd)
%
%   % Example #2 - Design a Chebyshev Type II parametric equalizer filter
%   % that cuts by 12 dB
%   d = fdesign.parameq('N,Flow,Fhigh,Gref,G0,GBW,Gst',...
%      4,.3,.5,0,-12,-10,-1);
%   Hd = design(d,'cheby2');
%   fvtool(Hd)
%
%   % Example #3 - Design a N=4 order audio lowpass (F0 = 0) shelving 
%   % filter with cutoff frequency of Fc = 0.25, quality factor Qa =10, and
%   % boost gain of G0 = 8 dB. 
%   d = fdesign.parameq('N,F0,Fc,Qa,G0',4,0,0.25,10,8);
%   Hd = design(d);
%   fvtool(Hd)
%
%   % Example #4 - Design a N=2 order audio parametric equalizer with
%   % center frequency F0 = 0.3, quality factor Qa = 5, reference gain of 0
%   % dB and boost gain of 10 dB.
%   d = fdesign.parameq('N,F0,Qa,Gref,G0',2,0.3,5,0,10);
%   Hd = design(d);
%   fvtool(Hd)
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN, FDESIGN/DESIGNOPTS.

%   Author(s): S. Orfanidis
%   Copyright 2008 The MathWorks, Inc.

this = fdesign.parameq;

set(this, 'Response', 'Parametric Equalizer');

this.setspecs(varargin{:});

capture(this);

% [EOF]
