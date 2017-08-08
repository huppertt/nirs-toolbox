function f = thisisparallel(Hd)
%THISISPARALLEL  True for filter with parallel sections.
%   THISISPARALLEL(Hd) returns 1 if filter Hd is composed of parallel sections,
%   and 0 otherwise. 
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2005 The MathWorks, Inc.

% This should be private

% The overloaded @parallel/isparallel will be called if Hd is paralleld.
f = false;

