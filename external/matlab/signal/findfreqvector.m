function [upn_or_w, upfactor, iswholerange, options, addpoint] = findfreqvector(varargin)
%FINDFREQVECTOR Define the frequency vector where the phase is evaluated.

%   Author(s): V.Pellissier
%   Copyright 2005 The MathWorks, Inc.


% Parse inputs
[n_or_w, options] = extract_norw(varargin{:});

% Add f=0 if not included in the frequency factor
addpoint = 0;
if length(n_or_w)>1,
    idx = find(n_or_w>=0);
    if isempty(idx),
        % Only negative frequencies
        addpoint = 1;
        n_or_w = [n_or_w,0];
    else
        idx = find(n_or_w<=0);
        if isempty(idx),
            % Only positive frequencies
            addpoint = -1;
            n_or_w = [0,n_or_w];
        end
    end
end

% Define the new N-point frequency vector where the frequency response is evaluated
[upn_or_w, upfactor, iswholerange] = getinterpfrequencies(n_or_w, varargin{:});
