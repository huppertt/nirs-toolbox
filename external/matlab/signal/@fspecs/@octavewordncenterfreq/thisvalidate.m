function [isvalid, errmsg, msgid] = thisvalidate(this)
%VALIDATE   Validate the specs

%   Copyright 2006-2012 The MathWorks, Inc.

isvalid = true;
errmsg  = '';
msgid   = '';

F0 = this.F0;
Fs = this.Fs;

if this.NormalizedFrequency,
     warning(message('signal:fspecs:octavewordncenterfreq:thisvalidate:NormalizedFrequency', 'NormalizedFrequency', 'normalizefreq'));
    Fs = 48000;
    F0 = F0*Fs/2;
end

validFreq = this.getvalidcenterfrequencies;
if isempty(find(F0 == validFreq, 1)),
    [~, idx] = min(abs(F0-validFreq));
    if this.NormalizedFrequency,
        F0 = validFreq(idx)/Fs*2;
        validfstr = num2str(validFreq/Fs*2);
    else
        F0 = validFreq(idx);
        validfstr = num2str(round(100*validFreq)/100);
    end
    this.F0 = F0;
    
    propName = 'F0';
    if isprop(this,'FromFilterDesigner') && this.FromFilterDesigner
      propName = signal.internal.filterdesigner.convertpropnames('F0');
    end    
    warning(message('signal:fspecs:octavewordncenterfreq:thisvalidate:InvalidF0',propName, validfstr,propName, num2str( F0 )));
end    

