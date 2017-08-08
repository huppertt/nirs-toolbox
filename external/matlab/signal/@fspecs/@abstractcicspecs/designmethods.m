function [d, isfull, type] = designmethods(this, varargin)
%DESIGNMETHODS   Return the design methods for this specification object.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

isfull = false;
type = 'fir';

d = {'multisection'};

% [EOF]
