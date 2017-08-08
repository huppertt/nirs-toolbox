function [streamsOnPool, useSubstreams, S, uuid] = unpackRNGscheme(RNGscheme)
%UNPACKRNGSCHEME extracts random stream info.
%
%   UNPACKRNGSCHEME is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.
%

%   Copyright 2010-2014 The MathWorks, Inc.

if ~isempty(RNGscheme)
    % RNGscheme fields
    uuid             = RNGscheme.uuid;
    useSubstreams    = RNGscheme.useSubstreams;
    streams          = RNGscheme.streams;
    useDefaultStream = RNGscheme.useDefaultStream;
    streamsOnPool    = RNGscheme.streamsOnPool;
    
    % Derived quantities
    if ~useDefaultStream && ~streamsOnPool
        % Serial computation with an explicitly supplied stream,
        % and/or if using Substreams.
        S = streams{1};
    else
        S = [];
    end
else
    % RNGscheme fields
    uuid = [];
    useSubstreams = false;
    streamsOnPool = false;
    
    % Derived quantities
    S = [];
end

