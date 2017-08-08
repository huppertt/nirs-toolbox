function Hflat = thisflatcascade(this,Hflat)
%THISFLATCASCADE Add singleton to the flat list of filters Hflat 

%   Copyright 2008 The MathWorks, Inc.

Hflat = [Hflat;this];