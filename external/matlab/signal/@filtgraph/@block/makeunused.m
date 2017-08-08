function [blk, delCoeffName] = makeunused(blk)
%MAKEUNUSED Sets a block as unused or DUMMY

%   Author(s): S Dhoorjaty
%   Copyright 1988-2004 The MathWorks, Inc.

setblocktype(blk,'DUMMY');
blk.mainParam = ''; blk.label = '';

% collect deleted coefficient name
delCoeffName = blk.coeffnames;