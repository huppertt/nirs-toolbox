function varargout = autoscale(this,x)
%AUTOSCALE   

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

if ~strcmpi(this.privArithmetic, 'fixed'),
    error(message('signal:dfilt:abstractfilter:autoscale:invalidArithmetic'));
end

% Verify that the structure support autoscale
verifyautoscalability(this);

w = warning('off');

if nargout>0,
    that = copy(this);
else
    that = this;
end

% Cache States 
states = that.States;

that.specifyall;

% Cache fipref
oldfipref = struct(fipref);
fipref('LoggingMode', 'on', 'DataTypeOverride', 'ScaledDoubles');
y = filter(that,x); 
R = qreport(that);

% Apply Safety margin
f = fieldnames(R);
for i=1:length(f),
    s.(f{i}).Min = R.(f{i}).Min;
    s.(f{i}).Max = R.(f{i}).Max;
end

% Autoscale
thisautoscale(that.filterquantizer,s,false);


% Restore fipref
fipref(oldfipref);

% Restore States
that.States = states;


if nargout>0,
    varargout{1} = that;
end

warning(w);


% [EOF]
