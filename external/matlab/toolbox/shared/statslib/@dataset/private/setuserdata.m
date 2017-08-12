function a = setuserdata(a,newdata)
%SETUSERDATA Set dataset array UserData property.

%   Copyright 2006-2011 The MathWorks, Inc.


if nargin < 2
    error(message('stats:dataset:setuserdata:TooFewInputs'));
end

a.props.UserData = newdata;
