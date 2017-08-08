function s = getdesignpanelstate(this, hfm)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

s = getdesignpanelstate(this.CurrentSpecs);

types = getfdatooltypes(this);

s.isDesigned = 1;
s.ResponseType = types{1};
s.SubType = types{2};
s.Tag = 'siggui.designpanel';
s.Version = 1;

% Design method must be specified by the FMETHOD object.
if nargin > 1
    sfm = getdesignpanelstate(hfm);
    s.DesignMethod = sfm.DesignMethod;
    if isfield(sfm, 'Components')
        s.Components = {s.Components{:}, sfm.Components{:}};
    end
end

% [EOF]
