function this = audioweighting(varargin)
%   AUDIOWEIGHTING  Construct an audio weighting filter designer.
%   D = FDESIGN.AUDIOWEIGHTING(SPECSTRING,VALUE1,VALUE2,...) constructs an audio
%   weighting filter designer D. Note that D is not the design itself. It only
%   contains the design specifications. In order to design the filter, invoke
%   the DESIGN method on D. For example (more examples below): 
%   D = fdesign.audioweighting('WT,Class','A',1); 
%   H = design(D); % H is a DFILT object
%
%   SPECSTRING is a string that determines what design specifications are used.
%   The full list of possible specifications is given below.
%
%   Different specifications may have different design methods available. Use
%   DESIGNMETHODS to get a list of design methods available for a given
%   SPECSTRING: 
%   designmethods(D)
%
%   VALUE1, VALUE2,... provide the values of the corresponding specifications.
%   Use get(D,'description') for a description of VALUE1, VALUE2, etc.
%
%   D = FDESIGN.AUDIOWEIGHTING(...,Fs) provides the sampling frequency in Hz 
%   used to design the filter. Fs must be a scalar trailing all other provided
%   values. Because audio weighting filter standards specify specific
%   attenuation values at a given frequency point, a filter designed for a
%   sampling frequency Fs only meets the specifications in the standard when
%   operating at that sampling frequency. If a sampling frequency is not
%   provided, the sampling frequency defaults to 48 kHz (80 KHz for the case of
%   ITU-R 468-4 designs). If you set normalized frequency to true using the
%   normalizefreq method, a warning is issued when the design method is invoked
%   and the default sampling frequency is used in the design.
%
%   The full list of possible values for SPECSTRING (not case sensitive) is:
%
%      'WT,Class' (default)
%      'WT'
%
%   where 
%       WT     - Weighting Type
%       Class  - Filter Class 
% 
%   The 'WT,Class' SPECSTRING will generate a filter designer for ANSI
%   S1.42-2001 weighting filters. Obtain a designer for an A-weighting filter by
%   setting the weighting type property to 'A'. Obtain a designer for a
%   C-weighting filter by setting the weighting type property to 'C'. The class
%   of these filters may be set to 1 or 2. The class value does not affect the
%   design. The class value is only used to provide a specification mask in
%   fvtool for the analysis of the filter design. Use the isspecmet method to
%   check that the filter design meets the specifications.
%
%   The 'WT' SPECSTRING generates a filter designer for a C-message (Bell System
%   Technical Reference 41009) weighting filter when you set the weighting type
%   property to 'Cmessage'. Obtain a filter designer for a Psophometer weighting
%   filter (ITU-T 0.41) by setting the weighting type property to 'ITUT041'.
%   Obtain an ITU-R 468-4 weighting filter designer by setting the weighting
%   type property to 'ITUR4684'.
%
%   D = FDESIGN.AUDIOWEIGHTING uses the default SPECSTRING ('WT,Class') and sets
%   the corresponding values to 'A', and 1. The default sampling frequency is 48
%   kHz.
%
%   % Example #1 - Design a class 2, A-weighting filter for a sampling 
%   % frequency of 44.1 KHz.
%   d  = fdesign.audioweighting('WT,Class','A',2,44.1e3)
%   Hd = design(d);
%   fvtool(Hd)
%
%   % Example #2 - Design a C-message IIR weighting filter for a sampling 
%   % frequency of 48 KHz.
%   d  = fdesign.audioweighting('WT','Cmessage',48e3)
%   Hd = design(d,'bell41009');
%   fvtool(Hd)
%
%   % Example #3 - Design a C-message FIR weighting filter for a sampling 
%   % frequency of 20 KHz.
%   d  = fdesign.audioweighting('WT','Cmessage',20e3)
%   Hd = design(d,'equiripple');
%   fvtool(Hd)
%
%   % Example #4 - Design a ITU-R 468-4 IIR weighting filter for a sampling 
%   % frequency of  44.1 KHz.
%   d  = fdesign.audioweighting('WT','ITUR4684',44.1e3)
%   Hd = design(d,'iirlpnorm');
%   fvtool(Hd)
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN, FDESIGN/DESIGNOPTS.

%   Copyright 2009 The MathWorks, Inc.

this = fdesign.audioweighting;

set(this, 'Response', 'Audio Weighting');

this.setspecs(varargin{:});

% [EOF]
