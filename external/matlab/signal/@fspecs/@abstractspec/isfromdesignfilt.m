function flag = isfromdesignfilt(this)
%ISFROMDESIGNFILT True if design comes from designfilt function

%   Copyright 2013 The MathWorks, Inc.

flag = isprop(this,'FromDesignfilt') && this.FromDesignfilt;