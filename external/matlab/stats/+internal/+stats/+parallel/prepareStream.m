function stream = prepareStream(iter,initialSubstream,S,usePool,useSubstreams,uuid)
%PREPARESTREAM readies the random stream object STREAM for a for/parfor iteration.
%
%   PREPARESTREAM gets the STREAM to use if the code is running on a parallel pool,
%   and there are multiple streams on different workers. It positions the 
%   Substream property if separate substreams are being used for each iterate.
%
%   PREPARESTREAM is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.

%   Copyright 2010-2014 The MathWorks, Inc.

if usePool
    S = internal.stats.parallel.workerGetValue(uuid);
end
if useSubstreams
    S.Substream = iter + initialSubstream - 1;
end
stream = S;
end %-prepareStream

