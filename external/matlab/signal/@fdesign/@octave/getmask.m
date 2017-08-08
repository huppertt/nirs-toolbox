function [F, A] = getmask(this, fcns, rcf, specs)
%GETMASK   Get the mask.

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

w = warning('off');
G = 10^(3/10);
b = this.BandsPerOctave;
omega_gt_1 = 1+((G^(1/(2*b))-1)/(G^(1/2)-1))*(G.^[0 1/8 1/4 3/8 1/2 1/2 1 2 3 4]-1);
omega_lt_1 = 1./omega_gt_1(end:-1:2);
FF = this.Currentspecs.F0*[omega_lt_1 omega_gt_1];

if strcmpi(this.Mask, 'class 0'),
    MMUp = [-75 -62 -42.5 -18 -2.3 .15 .15 .15 .15 .15 .15 .15 .15 .15 -2.3 -18 -42.5 -62 -75];
    MMLow = [NaN NaN NaN NaN -4.5 -4.5 -1.1 -.4 -.2 -.15 -.2 -.4 -1.1 -4.5 -4.5 NaN NaN NaN NaN];
elseif strcmpi(this.Mask, 'class 1'),
    MMUp = [-70 -61 -42 -17.5 -2 .3 .3 .3 .3 .3 .3 .3 .3 .3 -2 -17.5 -42 -61 -70];
    MMLow = [NaN NaN NaN NaN -5 -5 -1.3 -.6 -.4 -.3 -.4 -.6 -1.3 -5 -5 NaN NaN NaN NaN];
elseif strcmpi(this.Mask, 'class 2'),
    MMUp = [-60 -55 -41 -16.5 -1.6 .5 .5 .5 .5 .5 .5 .5 .5 .5 -1.6 -16.5 -41 -55 -60];
    MMLow = [NaN NaN NaN NaN -5.5 -5.5 -1.6 -.8 -.6 -.5 -.6 -.8 -1.6 -5.5 -5.5 NaN NaN NaN NaN];
else
    error(message('signal:fdesign:octave:getmask:InternalError'));
end
F = [FF NaN FF];
F = F*(fcns.getfs()/2);
if ~this.NormalizedFrequency
    F = F/(this.Fs/2);
end 
A = convertmagunits([MMUp NaN MMLow], 'db', 'linear', 'amplitude');
A = fcns.getarbmag(A);
warning(w);

% [EOF]
