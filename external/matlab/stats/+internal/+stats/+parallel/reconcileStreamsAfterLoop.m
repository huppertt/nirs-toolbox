function reconcileStreamsAfterLoop(RNGscheme,initialSubstreamOffset,niter)
%RECONCILESTREAMSAFTERLOOP updates random streams state after a parallelized loop.
%   RECONCILESTREAMSAFTERLOOP is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.

%   Copyright 2010-2014 The MathWorks, Inc.

    streams = RNGscheme.streams;
    if RNGscheme.useSubstreams
        streams{1}.Substream = initialSubstreamOffset + niter;
    else if RNGscheme.streamsOnPool
        poolStreams = internal.stats.parallel.retrieveFromPool(RNGscheme.uuid);
        for i=1:length(streams)
            streams{i}.State = poolStreams{i}.State;
        end
    end
end
