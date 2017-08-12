function a = scalarRepmat(a,m,n)
%SCACALRREPMAT Replicate a scalar dataset array

%   Copyright 2012 The MathWorks, Inc.


% Special handling for this case, because the repmat function
% uses linear indexing, which dataset does not allow.
a.nobs = m;
a.nvars = n;
a_data = a.data{1};
a_data = repmat(a_data,[m ones(1,ndims(a_data)-1)]);
a.data = repmat({a_data},1,n);

a.varnames = matlab.lang.makeUniqueStrings(repmat(a.varnames,1,n),{},namelengthmax);
if ~isempty(a.obsnames), a.obsnames = matlab.lang.makeUniqueStrings(repmat(a.obsnames,1,m),{},namelengthmax); end
if ~isempty(a.props.VarDescription), a.props.VarDescription = repmat(a.props.VarDescription,1,n); end
if ~isempty(a.props.Units), a.props.Units = repmat(a.props.Units,1,n); end
