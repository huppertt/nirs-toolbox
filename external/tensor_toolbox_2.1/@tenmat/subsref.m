function a = subsref(t,s)
%SUBSREF Subscripted reference for tenmat.
%
%   Examples
%   X(i,j) <-- returns the (i,j) entry in X
%   X.data <-- returns a 2D array of the data
%   X.tsize <-- returns the size original tensor
%   X.rdims <-- tensor dimensions that were mapped to rows
%   X.cdims <-- tensor dimensions that were mapped to columns 
%
%   See also TENMAT.
%
%MATLAB Tensor Toolbox.
%Copyright 2006, Sandia Corporation. 

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda. 
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2006) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: subsref.m,v 1.7 2006/08/21 21:04:40 bwbader Exp $

switch s(1).type
    case '.'
        switch s(1).subs
            case 'data'
                a = tt_subsubsref(t.data,s);
            case 'tsize'
                a = tt_subsubsref(t.tsize,s);
            case {'rindices','rdims'}
                a = tt_subsubsref(t.rindices,s);
            case {'cindices','cdims'}
                a = tt_subsubsref(t.cindices,s);
            otherwise
                error(['No such field: ', s.subs]);
        end
    case '()'
        a = t.data(s.subs{:});
    otherwise
        error('Invalid subsref into tenmat.')
end
end
