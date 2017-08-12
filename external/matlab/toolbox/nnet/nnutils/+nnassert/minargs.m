function minargs(numargs,minargs)
%NNASSERT_MINARGS Check number of arguments against minimum arguments.
%
%  NNASSERT_MINARGS(NARGIN,MINARGS) throws an error in the calling function
%  if NARGIN < MINARGS.

% Copyright 2010 The MathWorks, Inc.

if numargs < minargs
  throwAsCaller(MException(nnerr.tag('Args',2),'Not enough input arguments.'));
end
