function abstract_loadpublicinterface(this, s)
%ABSTRACT_LOADPUBLICINTERFACE   Load the public interface.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

base_loadpublicinterface(this, s);

% Make sure we force a copy of any FIs.
set(this, 'States', forcecopy(this, s.States));

% [EOF]
