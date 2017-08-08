function this = fdesign
%FDESIGN Digital Filter Design.
%   D = FDESIGN.<RESPONSE> returns a set of specifications that can be used
%   to design filters specified by RESPONSE, such as lowpass, highpass, and
%   many other types of filters. Type "help fdesign/responses" to get the
%   list of supported <a href="matlab:help fdesign/responses">responses</a>.
%
%   The design process creates the filter coefficients and associates a
%   particular filter structure to those coefficients. To design a filter
%   use the following command:
%   <a href="matlab:help fdesign/design">design</a>          -  Design the filter from the specifications.
%
%   Other functions that can aid the design process:
%   <a href="matlab:help fdesign/designmethods">designmethods</a>   -  Available design methods for the specifications given.
%   <a href="matlab:help fdesign/help">help</a>            -  Display help for a particular design method.
%   <a href="matlab:help fdesign/designopts">designopts</a>      -  Return a structure of design options to be used with <a href="matlab:help fdesign/design">design</a>.
%   <a href="matlab:help fdesign/designoptions">designoptions</a>   -  Show all design options available for a particular design.
%   <a href="matlab:help fdesign/validstructures">validstructures</a> -  Show valid filter structures for a particular design.
%
%   To process data through the filter, see <a href="matlab:help dfilt/filter">dfilt/filter</a> and <a href="matlab:help mfilt/filter">mfilt/filter</a>.
%
%   EXAMPLE - Design a minimum order lowpass filter with
%             % a normalized passband frequency of 0.2, 
%             % a stopband frequency of 0.22, 
%             % a passband ripple of 1 dB, 
%             % and a stopband attenuation of 60 dB. 
%             d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.2, 0.22, 1, 60);
%             f = design(d, 'equiripple');  % Design an FIR equiripple filter
%             info(f)                       % View information about filter
%
%             % Other designs can be performed for the same specifications
%             designmethods(d,'iir'); % List the available IIR design methods
%             f = design(d, 'ellip')  % Design an elliptic IIR filter (SOS)
%             fvtool(f)               % visualize various filter responses
%             input = randn(100,1);       
%             output = filter(f,input); % Process data through the elliptic filter.
%
%             % Many designs such as equiripple allow for various options to 
%             % be specified at design time. 
%
%   For more information enter "doc fdesign" at the MATLAB command line.
%
%   <a href="matlab:help signal">Signal Processing Toolbox TOC</a> 
%
%   See also FILTERBUILDER, FDESIGN/SETSPECS, FDESIGN/NORMALIZEFREQ.

%   Copyright 1999-2006 The MathWorks, Inc.

error(message('signal:fdesign:fdesign:Package', 'FDESIGN.*', 'h = fdesign.lowpass', 'help fdesign/responses'));


% [EOF]
