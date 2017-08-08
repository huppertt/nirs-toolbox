function y = sigcasttofloat(x, castType, fcnName, varName,datacheckflag)
%SIGCASTTOFLOAT Check if input x is floating point or numeric and then
%casts input to castType
% Inputs:
% x             - input data
% castType      - data type we want to cast to ('single','double')                 
% fcnName       - function name
% varName       - variable name
% datacheckflag - can be 'allowfloat' or 'allownumeric'. Default is
%                 'allowfloat'. When set to 'allowfloat' the function
%                 checks if data is floating point and then casts data to
%                 the specified castType type. When set to 'allownumeric'
%                 the function checks if data is numeric and then casts
%                 data to the specified castType type.
%
% Outputs:
% y            - cast output data

%   Copyright 2013 The MathWorks, Inc.

if nargin < 3
  fcnName = '';
  varName = '';
  datacheckflag = 'allowfloat';
end
if  nargin < 4
  varName = '';
  datacheckflag = 'allowfloat';
end
if nargin < 5
  datacheckflag = 'allowfloat';
end

signal.internal.sigcheckfloattype(x,'', fcnName, varName, datacheckflag);

y = cast(x,castType);
