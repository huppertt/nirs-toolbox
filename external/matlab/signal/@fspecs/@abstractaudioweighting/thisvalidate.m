function [isvalid, errmsg, msgid] = thisvalidate(this)
%THISVALIDATE   Validate the specs

%   Copyright 2009 The MathWorks, Inc.

isvalid = true;
errmsg  = '';
msgid   = '';

% If normalize frequency is true, warn and set the design Fs to the default Fs
% value. Otherwise, set the design Fs to the Fs specified by the user. The
% weighting standards specify attenuation values for specific frequency points
% in Hz. For this reason, the sampling frequency becomes a design parameter. 
if this.NormalizedFrequency,
     warning(message('signal:fspecs:abstractaudioweighting:thisvalidate:NormalizedFrequency', 'NormalizedFrequency', sprintf( '%g', this.DefaultFs )));
     
     this.ActualDesignFs = this.DefaultFs;          
else
    this.ActualDesignFs = this.Fs;
end

   
% [EOF]
