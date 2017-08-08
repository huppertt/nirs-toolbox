function [zall, pall, kall] = zplane(this, opts)
%ZPLANE Returns the zeroes and poles of filters

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin < 2,
    opts.showref  = false;
    opts.showpoly = false;
    opts.sosview  = [];
end

if length(this) > 1 || ~isa(this(1).Filter, 'dfilt.abstractsos')
    opts.sosview = [];
end

zall = {};
pall = {};
kall = {};

for indx = 1:length(this),
    cFilt = this(indx).Filter;
    if ispolyphase(cFilt) && opts.showpoly
        cFilt = polyphase(cFilt, 'objects');
    end

    if ~isempty(opts.sosview)
        cFilt = getfilters(opts.sosview, cFilt);
    end

    for jndx = 1:length(cFilt)
        [z, p, k] = zplane(cFilt(jndx));
        if iscell(z)
            z = [z{:}];
            z = z(:);
            p = [p{:}];
            p = p(:);
            end
            zall = {zall{:}, z};
            pall = {pall{:}, p};
            kall = {kall{:}, k};
        if opts.showref && isquantized(cFilt(jndx))
            [z, p, k] = zplane(reffilter(cFilt(jndx)));
            if iscell(z)
                z = [z{:}];
                z = z(:);
                p = [p{:}];
                p = p(:);
            end
            zall{end} = [nanpad(zall{end}, length(z)) z];
            pall{end} = [nanpad(pall{end}, length(p)) p];
            kall{end} = [kall{end} k];
        end
    end
end

% ------------------------------------------------------------------------
function p = nanpad(p, n)

p = [p(:); repmat(NaN, n-length(p), 1)];

% [EOF]
