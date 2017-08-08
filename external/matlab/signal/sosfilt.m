function X = sosfilt(SOS,X,dim)
%SOSFILT Second order (biquadratic) IIR filtering.
%   SOSFILT(SOS, X) filters the data in vector X with the second-order
%   section (SOS) filter described by the matrix SOS. The coefficients of
%   the SOS matrix must be expressed using an Lx6 second-order section
%   matrix where L is the number of second-order sections. If X is a
%   matrix, SOSFILT will filter along the columns of X. If X is a
%   multidimensional array, sosfilt operates on the first nonsingleton
%   dimension. SOS, X, or both can be double or single precision. If one
%   input is single precision then filtering is done with single precision
%   arithmetic.
%
%   SOSFILT uses a direct form II implementation to perform the filtering.
%
%   The SOS matrix should have the following form:
%
%   SOS = [ b01 b11 b21 a01 a11 a21
%           b02 b12 b22 a02 a12 a22
%           ...
%           b0L b1L b2L a0L a1L a2L ]
%
%   SOSFILT(SOS, X, DIM) operates along the dimension DIM.
%
%   % Example:
%   %   Filter a noisy signal with a lowpass IIR filter using sosfilt.
%   
%   Fs = 100;                                   % sample rate
%   t = 0:1/Fs:1;                               
%   x1 = sin(2*pi*t*3) + 0.1*randn(size(t));    % 3 Hz noisy sine wave
%   [z,p,k] = butter(4,5/(Fs/2));               % lowpass filter - 5 Hz cutoff
%   SOS = zp2sos(z,p,k);                        % second order sections matrix
%   y = sosfilt(SOS,x1);                        % filter the signal
%   plot(t,x1,'k',t,y,'g')
%   legend('Original Signal','Filtered Signal'); grid on
%
%   See also LATCFILT, FILTER, TF2SOS, SS2SOS, ZP2SOS, SOS2TF, SOS2SS, SOS2ZP.

%   Copyright 1988-2012 The MathWorks, Inc.
%

error(nargchk(2,3,nargin),'struct')

if nargin<3
  dim = [];
end

if ~isfloat(X)
  error(message('signal:sosfilt:VectMustBeFloat', 'X'));
end

if ~isfloat(SOS)
  error(message('signal:sosfilt:MatrixMustBeFloat', 'SOS'));
end

[m,n]=size(SOS);
if (m<1) || (n~=6),
  error(message('signal:sosfilt:InvalidDimensions'));
end

if isa(X,'single') || isa(SOS,'single')
  % Single precision filtering
  
  % SOSFILT uses a direct form II implementation to perform the filtering
  H = dfilt.df2sos(SOS);
  H.Arithmetic = 'single';
  X = filter(H,X,dim);
else
  % Double precision filtering
  
  h = SOS(:,[5 6 1:3]);
  for i=1:size(h,1),
    h(i,:)=h(i,:)./SOS(i,4);  % Normalize by a0
    h(i,[1 2]) = -h(i,[1 2]); % [-a1 -a2 b0 b1 b2]
  end
  h=h.';
  
  s = size(X);
  
  [X,perm,nshifts] = shiftdata(X,dim);
  s_shift = size(X); % New size
  X = reshape(X,size(X,1),[]); % Force into 2-D
  
  % Mex file will always filter along the columns
  X = sosfiltmex(h,X);
  
  % Convert back to the original shape
  X = reshape(X,s_shift); % Back to N-D array
  X = unshiftdata(X,perm,nshifts);
  X = reshape(X,s);
end