function t = subsasgn(t,s,b)
%SUBSASGN Subscripted assignement for ktensor.
%
%   See also KTENSOR.
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
% $Id: subsasgn.m,v 1.5 2006/08/29 00:39:09 tgkolda Exp $

switch s(1).type
    case '.'
        switch s(1).subs
            case 'lambda'
                if length(s) == 1
                    t = ktensor(b, t.u);
                else
                    newlambda = subsasgn(t.lambda, s(2:end), b);
                    t = ktensor(newlambda, t.u);
                end
            case {'u','U'}
                if length(s) == 1
                    t = ktensor(t.lambda, b);
                else
                    tmpu = subsasgn(t.u, s(2:end), b);
                    t = ktensor(t.lambda, tmpu);
                end
            otherwise
                error(['No such field: ', s(1).subs]);
        end
    case '()'
        error('Cannot change individual entries in a ktensor.')
    case '{}'
        new_s(1).type = '.';
        new_s(1).subs = 'u';
        new_s(2:length(s)+1) = s;
        t = subsasgn(t, new_s, b);
    otherwise
        error('Invalid subsasgn.');
end


