function [params, values, descs, args] = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

[params, values, descs, args] = firceqrip_genmcode(h, d);

% [EOF]
