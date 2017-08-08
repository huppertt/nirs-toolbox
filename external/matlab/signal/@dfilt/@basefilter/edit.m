function edit(this)
%EDIT   Edit the dfilt.

%   Author(s): J. Schickler
%   Copyright 2005-2006 The MathWorks, Inc.

hfdesign = getfdesign(this);
if isempty(hfdesign)
    
    % If there is no associated FDESIGN object, just bring up the
    % coefficient editor.
    h = FilterDesignDialog.CoeffEditor(this);
else
    % Construct the dialog given the fdesign object.
    h = feval(getdialogconstructor(hfdesign));

    vname = inputname(1);

    setGUI(h, this)
end


hDlg = DAStudio.Dialog(h);

l = handle.listener(h, 'DialogApplied', ...
    @(hObj, ed) assignin('base', vname, design(h)));
schema.prop(hDlg, 'DialogAppliedListener', 'handle.listener');
set(hDlg, 'DialogAppliedListener', l);

% [EOF]
