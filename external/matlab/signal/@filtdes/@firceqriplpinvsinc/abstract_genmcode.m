function [params, values, descs, args] = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[params, values, descs, args] = firceqrip_genmcode(h, d);

params = {params{:}, 'isincffactor', 'isincpower'};
values = {values{:}, getmcode(d, 'invsincfreqfactor'), getmcode(d, 'invsincpower')};
descs  = {descs{:}, 'Inverse Sinc Frequency Factor', 'Inverse Sinc Power'};

args = sprintf('%s, ''invsinc'', [isincffactor isincpower]', args);

% [EOF]
