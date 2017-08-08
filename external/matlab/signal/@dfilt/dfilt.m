function varargout=dfilt(varargin) %#ok<STOUT>
%DFILT  Digital Filter Implementation.
%   Hd = DFILT.<STRUCTURE>(COEFF,...) returns a discrete-time filter Hd that
%   associates the coefficients COEFFS with a particular filter STRUCTURE.
%   Coefficients are specified with one or more inputs depending on the
%   structure. 
%   Type "help dfilt/structures" to get the complete list of <a href="matlab:help dfilt/structures">structures</a>.
%
%   Digital filters are equipped with a variety of functions in categories of 
%   <a href="matlab:help dfilt/freqtransform">frequency transformation</a>, <a href="matlab:help dfilt/analysis">analysis</a>, <a href="matlab:help dfilt/simulation">simulation</a> , <a href="matlab:help dfilt.basefilter.fixedpointoptimizations">fixed-point optimizations</a>
%   and <a href="matlab:help dfilt/codegeneration">code generation</a>. 
%   (Type help dfilt/<function category> for more details, e.g., help dfilt/analysis)
%
%   The functions most commonly used with digital filters are:
%   <a href="matlab:help dfilt/filter">filter</a>        - Execute ("run") the discrete-time filter.
%   <a href="matlab:help dfilt/freqz">freqz</a>         - Compute the frequency response of the discrete-time filter.
%   <a href="matlab:help dfilt/realizemdl">realizemdl</a>    - Generate a Simulink subsystem.   (Simulink Required) 
%   <a href="matlab:help dfilt/block">block</a>         - Generate a Digital Filter block. (Simulink and DSP System Toolbox Required)
%   <a href="matlab:help dfilt/generatehdl">generatehdl</a>   - Generate HDL code.               (Filter Design HDL Coder Required)
%   <a href="matlab:help dfilt/sysobj">sysobj</a>        - Generate a filter System object. (DSP System Toolbox Required)
%
%   Notice that the DSP System Toolbox, along with the Fixed-Point
%   Designer, enables single precision floating-point and fixed-point
%   support for most DFILT structures.
%
%   % EXAMPLE: Design and construct a direct-form FIR lowpass filter and
%   % analyze its various responses
%   b = firls(80,[0 .4 .5 1],[1 1 0 0],[1 10]);
%   Hd = dfilt.dffir(b)
%   fvtool(Hd) % Analyze filter
%
%   For more information enter "doc dfilt" at the MATLAB command line.
%
%   <a href="matlab:help signal">Signal Processing Toolbox TOC</a> 
%
%   See also FDESIGN. 

%   Copyright 1988-2013 The MathWorks, Inc.

error(message('signal:dfilt:dfilt:Package'));
