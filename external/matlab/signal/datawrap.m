function [x,msg] = datawrap(x,nfft)
%DATAWRAP Wrap input data modulo nfft.
%   DATAWRAP(X,NFFT) wraps the vector X modulo NFFT.
%   
%   The operation consists of dividing the vector X into segments each of
%   length NFFT (possibly padding with zeros the last segment).  Subsequently,
%   the length NFFT segments are added together to obtain a wrapped version of X.

%   Author(s): R. Losada 
%   Copyright 1988-2004 The MathWorks, Inc.

nx = size(x, 2);
msg = '';
if all(size(x)>1),
   error(message('signal:datawrap:InvalidInput'));
end

% Reshape into multiple columns (data segments) of length nfft.
% If insufficient data points are available, zeros are appended.
% Sum across the columns (data segments).
x = sum(buffer(x,nfft),2);
% Reshape vector as necessary:
if (nx~=1), x=x.'; end

% [EOF] datawrap.m
