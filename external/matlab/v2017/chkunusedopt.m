function chkunusedopt(arglist)
%CHKUNUSEDOPTION error if any string options are found in the arguments
%   This function is for internal purposes only and may be removed in a
%   future release.
%
%   See also GETMUTEXCLOPT.

%   Copyright 2015 The MathWorks, Inc.

idx = find(cellfun(@ischar,arglist));
if any(idx)
  error(message('signal:chkunusedopt:UnrecognizedOption',arglist{idx(1)}));
end