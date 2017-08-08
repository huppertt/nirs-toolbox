function this = octave(varargin)
%OCTAVE   Construct an octave filter designer.
%   D = FDESIGN.OCTAVE(L) constructs a fractional-octave-band filter
%   designer D with L bands per octave. L defaults to 1.
%
%   D = FDESIGN.OCTAVE(L,MASK) initializes the filter designer 'Mask'
%   property to MASK.  MASK is one of the following strings: 'Class 0',
%   'Class 1', 'Class 2'. Notice that MASK is not a design parameter. It is
%   just used to draw the specification mask in FVTool.
%
%   D = FDESIGN.OCTAVE(L,MASK,SPEC) initializes the filter designer
%   'Specification' property to SPEC.  SPEC is one of the following strings
%   and is not case sensitive:
%
%       'N,F0'  
%
%  where 
%       N     - Filter Order
%       F0    - Center Frequency
%
%   Use GET(D, 'DESCRIPTION') for a description of SPEC. By default, all
%   frequency specifications are assumed to be in normalized frequency
%   units. 
%
%   D = FDESIGN.OCTAVE(L, MASK, N, F0) uses the  default SPEC ('N,F0') and
%   sets the filter order and center frequency.
%
%   D = FDESIGN.OCTAVE(...,Fs) specifies the sampling frequency (in Hz). In
%   this case, the center frequency is also in Hz and must be between 20Hz
%   and 20kHz (audio range).
%
%   % Example #1 - Design an 6th order, octave-band class 0 filter with:
%   %              a center frequency of 1000 Hz and,
%   %              a sampling frequency of 44.1kHz.
%   d = fdesign.octave(1,'Class 0','N,F0',6,1000,44100)
%   Hd = design(d)
%   fvtool(Hd,'Fs',44100,'FrequencyScale','log')
%
%   % Example #2 - Design a 6th order, 1/3-octave-band Class 2 filter with:
%   %              a center frequency of 10 KHz,
%   %              a sampling frequency of 48 kHz,
%   %              a Direct-Form I, Second-Order Sections structure.
%   d = fdesign.octave(3,'Class 2','N,F0',6,10000,48000)
%   design(d,'FilterStructure','df1sos')
%
%   For more information about octave-band filters, see the
%   <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\octavedemo.html'])">Octave-Band and Fractional-Octave-Band Filters</a> demo. 
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN.

%   Author(s): V. Pellissier
%   Copyright 2006-2010 The MathWorks, Inc.

this = fdesign.octave;
set(this, 'Response', 'Octave and Fractional Octave');
if nargin>0,
    this.BandsPerOctave = varargin{1};
    varargin(1) = [];
end
if nargin>1,
    this.Mask = varargin{1};
    varargin(1) = [];
end
    
setspecs(this, varargin{:});

% [EOF]
