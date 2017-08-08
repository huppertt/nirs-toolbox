function [upn_or_w, upfactor, iswholerange, do_transpose] = getinterpfrequencies(n_or_w, varargin)
%GETINTERPFREQUENCIES  Define the interpolation factor for phasez and zerophase.
%   [UPN_OR_W, UPFACTOR, ISWHOLERANGE] = GETINTERPFREQUENCIES(N_OR_W, VARARGIN) returns
%   the nfft (respectively w frequencies vector) UPN_OR_W to pass to freqz that is 
%   greater than a threshold (2^13), the upsampling factor UPFACTOR and a 
%   the ISWHOLERANGE boolean. 

%   Author(s): V.Pellissier, R. Losada
%   Copyright 1988-2009 The MathWorks, Inc.

% Minimum number of point where the frequenciy response will be evaluated.
threshold = 2^13;

% Determine if the whole range is needed
iswholerange = 0;
if nargin>2 && any(strcmpi('whole', varargin)),
    iswholerange = 1;
end

isn = 0;
N = length(n_or_w);
if length(n_or_w)==1,
    isn = 1;
    N = n_or_w;
end
 
% Default values
upfactor = 1;

% Compute the upfactor
if N<threshold,
    upfactor = ceil(threshold/N);
    if iswholerange,
        upfactor = 2*upfactor;
    end
else
    if iswholerange,
        upfactor = 2;
    end
end
    
do_transpose = false;
if isn,
    upn_or_w = N*upfactor;
else
    % Interpolate w if needed
    w = n_or_w(:);
    
    if upfactor == 1
        upn_or_w = w;
    else
        
        % Originally using interp, but that required 2*L+1 frequencies passed
        % in. This was broken for case when a two element frequencies vector
        % was passed in.
        %         upn_or_w = interp(w, upfactor);
        
        % Add one sample to the end of the w vector so that we can include
        % the last frequency specified by the user in the interpolated w
        % vector.
        w(end+1) = 2*w(end) - w(end-1);
        n = length(w);
        
        % Preallocate for upn_or_w vector
        upn_or_w = zeros(upfactor*(n-1),1);
        for idx = 0:n-2,
            beginIDX = (idx*upfactor)+1;
            endIDX = (idx+1)*upfactor;
            upn_or_w(beginIDX:endIDX,1) = linspace(w(idx+1),w(idx+2), upfactor);
        end
        if size(n_or_w,2)==1
            do_transpose = true;
        end
    end
end

% [EOF]
