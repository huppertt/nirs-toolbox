function Hd = fftfir(num,L)
%FFTFIR Overlap-add FIR.
%   Hd = DFILT.FFTFIR(NUM,L) constructs a discrete-time FIR filter object
%   for filtering using the overlap-add method. 
%
%   NUM is a vector of numerator coefficients.
%
%   L is the length of each block of input data used in the filtering.
%
%   The number of FFT points is given by L+length(NUM)-1. It may be
%   advantageous to choose L such that the number of FFT points is a power
%   of two.
%  
%   Note that one usually does not construct DFILT filters explicitly.
%   Instead, one obtains these filters as a result from a design using <a
%   href="matlab:help fdesign">FDESIGN</a>. 
%
%   % EXAMPLE #1: Direct instantiation
%   b = [0.05 0.9 0.05];
%   len = 50;
%   Hd = dfilt.fftfir(b,len)
%
%   % EXAMPLE #2: Design an equiripple lowpass filter with default specifications
%   Hd = design(fdesign.lowpass, 'equiripple', 'Filterstructure', 'fftfir');
%   fvtool(Hd)        % Analyze filter
%   x = randn(100,1); % Input signal
%   y = filter(Hd,x); % Apply filter to input signal
%
%   See also DFILT/STRUCTURES
  
%   Copyright 1988-2012 The MathWorks, Inc.
Hd = dfilt.fftfir;

Hd.FilterStructure = 'Overlap-Add FIR';

if nargin>=1
    Hd.Numerator = num;
end

if nargin < 2,
    % Set a default blocklength
    % Don't use factoryValue so overload set runs
    L = 100;
end

Hd.BlockLength = L;
