function p = get(a,name)
%GET Get a dataset array property.
%   GET(A) prints the list of properties for the dataset A.
%
%   S = GET(A) returns all properties of A in a scalar structure.
%
%   P = GET(A,'PropertyName') returns the property 'PropertyName'.
%
%   See also DATASET/SET, DATASET/SUMMARY.

%   Copyright 2006 The MathWorks, Inc. 


if nargin == 1
    s = a.props; s.ObsNames = a.obsnames; s.VarNames = a.varnames;
    if nargout == 1
        p = s;
    else
        disp(s);
    end
    
elseif nargin == 2
    if iscellstr(name)
        p = cell(1,numel(name));
        for i = 1:length(p)
            p{i} = getproperty(a,name{i});
        end
    else
        p = getproperty(a,name);
    end
end
