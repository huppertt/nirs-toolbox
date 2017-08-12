function flag = sigcheckfloattype(x, dataType, fcnName, varName, datacheckflag)
%SIGCHECKFLOATTYPE Check if input x is floating point or numeric and of the
%expected data type
%
% Inputs:
% x             - input data
% dataType      - data type we want to check ('single','double','int8',...)
%                 if set to empty ('') then we do not check for a specific
%                 data type. We only check if data is floating point or
%                 numeric depending on the datacheckflag input.
% fcnName       - function name
% varName       - variable name
% datacheckflag - can be 'allowfloat' or 'allownumeric'. Default is
%                 'allowfloat'. When set to 'allowfloat' the function
%                 checks if data is floating point and then checks if data
%                 is of the specified dataType type. When set to
%                 'allownumeric' the function checks if data is numeric and
%                 then checks if data is of the specified dataType type.
%
% Outputs:
% flag          - true if data is of type dataType


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

if strcmpi(datacheckflag,'allowfloat')
  typeCheck = isfloat(x);
  expType = 'double/single';
elseif strcmpi(datacheckflag,'allownumeric')
  typeCheck = isnumeric(x);
  expType = 'numeric';
else
  error(message('signal:sigcheckfloattype:InvalidDataCheckFlag'));    
end

if ~typeCheck
  if ~isempty(fcnName)
    if ~isempty(varName)
      error(message('signal:sigcheckfloattype:InvalidInput',...
        varName, fcnName, expType, class(x)));    
    else
      error(message('signal:sigcheckfloattype:InvalidInput1',...
        fcnName, expType, class(x)));    
    end
  else
    if ~isempty(varName)
      error(message('signal:sigcheckfloattype:InvalidInput2',...
        varName, expType, class(x)));
    else
      error(message('signal:sigcheckfloattype:InvalidInput3',...
        expType, class(x)));
    end
  end    
end

flag = isa(x,dataType);
