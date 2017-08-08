function [n_or_w, options, do_transpose] = extract_norw(varargin)
%EXTRACT_NORW Deal optional inputs of phasez and zerophase.
%   [N_OR_W, OPTIONS]=EXTRACT_NORW(VARARGIN) return the third input of freqz N_OR_W
%   that can be either nfft or a vector of frequencies where the frequency response 
%   will be evaluated.

%   Author(s): V.Pellissier, R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

% Default values
n_or_w = 512;
options = {};
do_transpose = false;

if nargin > 0 & isnumeric(varargin{1}) & isreal(varargin{1}),
    if ~isempty(varargin{1})
        n_or_w = varargin{1};
        if length(n_or_w)>1 && size(n_or_w,1)==1 
           do_transpose = true;
        end
        % Force row vector
        n_or_w = n_or_w(:).';
    end
    if nargin>1,
        options = {varargin{2:end}};
    end
else
    options = {varargin{:}};
end


% [EOF]
