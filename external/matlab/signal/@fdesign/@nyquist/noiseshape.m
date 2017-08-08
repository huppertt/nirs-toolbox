function Hns = noiseshape(this,Hd,WL,args)
%NOISESHAPE Noise-shape the FIR filter Hd

% This should be a private method

%   Copyright 2008 The MathWorks, Inc.

Hns = reffilter(Hd);
NBands = this.Band;
b = Hns.Numerator;
nb = length(b);
p = firpolyphase(b,NBands);

criticalband = [0 NBands*args.Fpass];
F = criticalband;
A = [1/NBands 1/NBands];

% Note for the halfband case, we decrease the wordlength by one in
% order to later accommodate for the 1/L coefficient which is not
% being used for the halfband noise-shaping optimization
WL = WL-1;

npoly = size(p,2);
% Find pure delay phase index
puredelay = find(sum(p==0,2)==npoly-1);
for i=1:NBands,
    if i~=puredelay,
        % Call super method
        linearphase = signalpolyutils('islinphase',p(i,:),1);
        nsres = supernoiseshape(this,p(i,:),linearphase,WL,criticalband,F,A,args);
        p(i,:) = nsres.filters.bns;
    end
end

% Maintain phase linearity
b = p(1:nb);
if islinphase(Hd),
    switch firtype(Hd),
        case 1,
            m = (nb-1)/2;
            b = [b(1:m), b(m+1), b(m:-1:1)];
        case 2,
            m = nb/2;
            b = [b(1:m), b(m:-1:1)];
        case 3,
            m = (nb-1)/2;
            b = [b(1:m), 0, -b(m:-1:1)];
        case 4,
            m = nb/2;
            b = [b(1:m), -b(m:-1:1)];
    end
end
Hns.Numerator = b;

% [EOF]