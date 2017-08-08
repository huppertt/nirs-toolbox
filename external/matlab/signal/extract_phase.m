function [phi,w] = extract_phase(h,upn_or_w,iswholerange,upfactor,w,addpoint)
%EXTRACT_PHASE Extract phase from frequency response  

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

% When h==0, the phase is not defined (introducing NaN's)
h = modify_fresp(h);

% Unwrap the phase
phi = unwrap_phase(h,upn_or_w,iswholerange);

% Downsample
phi = downsample(phi, upfactor);
w = downsample(w, upfactor);

% Remove additional point
if addpoint==1,
    phi(end)=[];
    w(end)=[];
elseif addpoint==-1,
    phi(1)=[];
    w(1)=[];
end

%-------------------------------------------------------------------------------
function h = modify_fresp(h)
% When h==0, the phase is not defined (introducing NaN's)

tol = eps^(2/3);
ind = find(abs(h)<=tol);
if ~isempty(ind);
    h(ind)=NaN;
end

%-------------------------------------------------------------------------------
function phi = unwrap_phase(h,w,iswholerange)

idx = find(w<0);
if isempty(idx),
    % Range only positive frequencies
    phi=unwrap(angle(h));
else
    idx = idx(end);
    % Unwrap negative frequencies
    phi_n=unwrap(angle(h(idx:-1:1)));
    if idx<length(w),
        % Unwrap positive frequencies
        phi_p=unwrap(angle(h(idx+1:end)));
    else
        phi_p = [];
    end
    phi=[phi_n(end:-1:1);phi_p];
end

% [EOF]
