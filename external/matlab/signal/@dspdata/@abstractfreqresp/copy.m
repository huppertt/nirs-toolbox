function h = copy(this)
%COPY   Copy the object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% Need to set all these properties at once since they're related.
propname = getrangepropname(this);
proplist = {...
    this.Data, ...
    this.Frequencies,...
    propname,get(this,propname),...
    'CenterDC',this.centerDC};

if this.NormalizedFrequency,
    h = feval(this.class,proplist{:});
    h.privFs  = getfs(this); % Store Fs in case it's not the default value.
else
    h = feval(this.class,proplist{:},'Fs',getfs(this));
end

h.Metadata                 = copy(this.Metadata);
h.privNormalizedFrequency  = this.NormalizedFrequency;


% [EOF]
