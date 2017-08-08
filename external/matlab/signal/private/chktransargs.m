function [x, t, n, d] = chktransargs(needDelay, x, varargin)
%CHKTRANSARGS check arguments for bilevel waveform transitions
%
%   Validates and extracts (numeric) inputs of the form:
%       (X,     <optional string arguments>)
%       (X, FS, <optional string arguments>)
%       (X, T,  <optional string arguments>)
%
%   If the needDelay argument is true then extract inputs of the form:
%       (X, D,     <optional string arguments>)
%       (X, FS, D, <optional string arguments>)
%       (X, T, D,  <optional string arguments>)
%
%   It will return the sample values, x, and sample instants, t,
%   that correspond to the input vector as well as the index, n, into
%   varargin that specifies the first [unhandled] string argument.
%   It will also return the delay parameter, d, if it is requested.
%
%   This function is for internal use only. It may be removed in the future.
   
%   Copyright 2011 The MathWorks, Inc.

  try
    validateattributes(x,{'double'},{'real','finite','vector'}, '','X');
    if (numel(x)<2)
      error(message('signal:chktransargs:MustBeMultiElementVector','X'));
    end
    x = x(:);
    
    n = getNumExtraArgs(varargin{:});
    
    if n>1+needDelay
      error(message('signal:chktransargs:TooManyNumericArgs',2+needDelay));
    end
    
    if needDelay
      if n == 0
        error(message('signal:chktransargs:MissingParameter','D'));
      end
      d = getDelayParameter(varargin{n});
    end
    t = getTimeVector(x, varargin{1:n-needDelay});
    n = n+1;
  catch ex
    throwAsCaller(ex);
  end
end

function n = getNumExtraArgs(varargin)
  % get the number of extra numeric arguments for transition measurements.
  n = 0;
  while n<nargin && isnumeric(varargin{n+1})
    n = n + 1;
  end
end

function t = getTimeVector(x, varargin)
  if nargin==1
    t = (1:numel(x))';
  elseif isscalar(varargin{1})
    fs = varargin{1};
    validateattributes(fs,{'double'},{'real','finite','positive'},'','FS');
    t = (0:numel(x)-1)' / fs;
  else
    t = varargin{1};
    validateattributes(t,{'double'},{'real','finite','vector'},'','T');
    t = t(:);
    if numel(x) ~= numel(t)
      error(message('signal:chktransargs:LengthMismatch'));
    elseif ~issorted(t) && ~isempty(diff(t)==0)
      error(message('signal:chktransargs:MustBeStrictlyIncreasing', 'T'));
    end
  end
end

function d = getDelayParameter(d)
  validateattributes(d,{'double'},{'real','finite','positive','scalar'},'','D');
end
