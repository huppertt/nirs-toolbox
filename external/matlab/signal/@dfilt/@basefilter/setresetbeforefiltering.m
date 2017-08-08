function flag =setresetbeforefiltering(Hd, flag)
%SETRESETBEFOREFILTERING Set function of the ResetBeforeFiltering property.

%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

Hd.PersistentMemory = ~strcmpi(flag,'on');

