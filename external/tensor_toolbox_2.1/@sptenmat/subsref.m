function a = subsref(t,s)
%SUBSREF Subscripted reference for a sptenmat.
%
%   Examples
%   A.subs <-- returns the nonzero values as an array
%   A.vals <-- returns the corresponding 2D subscripts
%   A.tsize <-- returns the size original tensor
%   A.rdims <-- tensor dimensions that were mapped to rows
%   A.cdims <-- tensor dimensions that were mapped to columns 
%
%   See also SPTENMAT.
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
% $Id: subsref.m,v 1.5 2006/08/21 21:04:39 bwbader Exp $

switch s(1).type    
    case '.'
        switch s(1).subs
            case 'vals'
                a = tt_subsubsref(t.vals,s);
	    case 'tsize'
		a = t.tsize;
	    case 'rdims'
		a = t.rdims;
	    case 'cdims'
		a = t.cdims;
	    case 'subs' 
                a = tt_subsubsref(t.subs,s);
            otherwise
                error(['No such field: ', s.subs]);
        end
    otherwise
        error('Invalid subsref into tenmat.')
end
