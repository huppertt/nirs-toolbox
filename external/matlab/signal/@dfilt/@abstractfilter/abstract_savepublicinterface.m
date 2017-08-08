function s = abstract_savepublicinterface(this)
%ABSTRACT_SAVEPUBLICINTERFACE   Save the public interface.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

s = base_savepublicinterface(this);

s.States = this.States;

% [EOF]
