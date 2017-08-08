function [resp, w] = completefreqresp(resp, fs, wmin, wmax, opts)
%COMPLETEFREQRESP Complete the frequency response for the specified range
%                 This is an undocumented function.

%   Inputs:
%       resp    -   A filter response from 0 to 2*pi
%       fs      -   The sampling Frequency of that response
%       fmin    -   The minimum frequency to return
%       fmax    -   The maximum frequency to return
%       opts    -   Structure with 3 fields:
%                    -periodicity 
%                    -flip 
%                    -shift 
%
%   Example #1:
%
%   [b,a] = butter(5,.5);
%   G = dfilt.df2t(b,a);
%   h = freqz(G, 512, 'whole');
%   [h, w] = completefreqresp(h, 90, -150, 150);
%   plot(w, 20*log10(abs(h)))
%
%   Example #2:
%
%   [b,a] = butter(5,.5);
%   G = dfilt.df2t(b,a);
%   p = phasez(G, 512, 'whole');
%   [pu, w] = completefreqresp(p, 90, -150, 150);
%   opts.periodicity = 2;
%   [ps, w] = completefreqresp(p, 90, -150, 150, opts);
%   plot(w, pu, w, ps)

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

error(nargchk(4,5,nargin,'struct'));

if nargin < 5
    opts = struct('periodicity', 2, 'flip', 0, 'shift', 0);
else
    if ~isfield(opts, 'periodicity'), opts.periodicity = 2; end
    if ~isfield(opts, 'flip'),        opts.flip        = 0; end
    if ~isfield(opts, 'shift'),       opts.shift       = 0; end
end

lorig  = length(resp);
period = opts.periodicity;

% Make sure the response is a column vector, this makes calculations easier.
if size(resp,1) == 1 || size(resp, 2) == 1
    resp = resp(:);
end
fsmax    = max(abs([wmin wmax]));
numresps = ceil(fsmax/(fs*period/2));

if opts.periodicity == 4,
    if opts.flip, m = -1;
    else          m = 1; end
    resp = [resp; flipud(resp)*m];
end

% If the minimum frequency is less than 0 we need to have information
% before zero.
if wmin < 0
    resp = repmat(resp, 2*numresps, 1);

    w = linspace(0, fs*(numresps*period/2-1/lorig), length(resp)/2)';
    w = [w-fs*numresps*period/2; w];
else
    resp = repmat(resp, numresps, 1);

    w = linspace(0, fs*(numresps*period/2-1/lorig), length(resp))';
end

% If there is a shift, apply it
if opts.shift ~= 0,
    
    % Build a shift value for each of the responses calculated above
    if wmin < 0
        shift = opts.shift*numresps:-opts.shift:-opts.shift*(numresps-1);
        numresps = numresps*2;
    else
        shift = 0:-opts.shift:-opts.shift*(numresps-1);
    end
    
    % Expand the shift vector for the length of the response (it is now a
    % matrix)
    shift = repmat(shift, length(resp)/numresps, 1);
    
    % Reshape the shift matrix so that it is a vector of the same length as
    % the response.
    shift = reshape(shift, length(resp), 1);

    % Add the shift to the response.
    resp  = repmat(shift, 1, size(resp, 2)) + resp;
end

% Exclude the highest frequency asked for
lindx       = find(w >= wmax);
w(lindx)    = [];
resp(lindx, :) = [];

% Include the lowest frequency asked for
findx       = find(w < wmin);
w(findx)    = [];
resp(findx, :) = [];

% [EOF]
