function [esttype, arglist] = psdesttype(validtypes,defaulttype,arglist)
%PSDESTTYPE - return the PSD estimation type option and remove it from the
%argument list
%
%   validtypes  - a cell array of valid estimator types
%                 (e.g. {'power','ms','psd'})
%
%   defaulttype - the default type to use if no type is found
%
%   varargin    - the input argument list
%
%   Errors out if different estimation types are specified in varargin.    

%   Copyright 2012-2013 The MathWorks, Inc.

esttype = defaulttype;
found = false;

for i=1:numel(validtypes)
  matches = find(strcmpi(validtypes{i},arglist));
  if ~isempty(matches)
    if ~found
      found = true;
      esttype = validtypes{i};
      arglist(matches) = [];
    else
      error(message('signal:psdoptions:ConflictingEstTypes', ...
                    esttype,validtypes{i}));
    end
  end
end
    
    
    