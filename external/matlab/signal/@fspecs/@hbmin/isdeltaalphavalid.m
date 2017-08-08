function alphaflag = isdeltaalphavalid(this,originalalpha,newalpha)
%ISDELTAALPHAVALID True if the object is deltaalphavalid

%   Copyright 2009 The MathWorks, Inc.

% We'd like to control how much alpha can vary according to the stop band
% attenuation.  We don't want to keep trying alpha if the spec just cannot
% be met.  We don't want to change alpha any more once its effect on
% stopband attenuation reach a certain point.

% Reference
% [1] Oppenheim, Discrete-Time Signal Processing, 1989, page 453

% stop trying alpha if its effect on Astop already reaches 3dB limit 
deltaAstop = 3; 
deltaalpha = newalpha - originalalpha;

Astop = this.Astop;
if Astop <= 21
    % in this case, alpha is 0, so we won't be able to change 
    alphaflag = false;
elseif (21 < Astop) && (Astop <= 50)
    % alpha = 0.5842(Astop-21)^0.4+0.07886(Astop-21)
    %       < sqrt(Astop-21)
    % Astop > alpha^2 + 21
    % deltaAstop > 2*alpha*deltaalpha
    alphaflag = deltaAstop > 2*deltaalpha*originalalpha;
else
    % alpha = 0.1102(Astop-8.7)
    % Astop > 9*alpha + 8.7
    % deltaAstop > 9*deltaalpha 
    alphaflag = deltaAstop > 9*deltaalpha;
end


% [EOF]
