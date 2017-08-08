function valid = iscompatibleRNGscheme(useParallel,RNGscheme)
%ISCOMPATIBLERNGSCHEME
%   ISCOMPATIBLERNGSCHEME is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.
%
%   ISCOMPATIBLERNGSCHEME tests 
%   (1) whether an RNGscheme is valid with the current parallel lpool state
%   (2) whether the RNGscheme is valid with the parallel computation flag.
%
%   In the case of (1) the RNGscheme is marked invalid if 
%   RNGscheme.streamsOnPool == true AND RNGscheme.streams is not the
%   same length as the parallel pool size.  In this case, we can determine 
%   that the RNGscheme is forever invalid because the parallel pool for which 
%   it was constructed has been closed (albeit perhaps replaced with a 
%   parallel pool of different size).  This test will NOT detect the case where
%   a parallel pool was reopened of the same size: that failure will occur 
%   when executable code tries to retrieve a value on a worker using 
%   RNGscheme.uuid as a retrieval key.
%
%   In the case of (2) the RNGscheme is invalid if
%   RNGscheme.streamsOnPool is true (in which case parallel computation is 
%   required), but useParallel is false.  It is also invalid if 
%   useParallel is true but the contents of RNGscheme don't support parallel 
%   computation (ie, ~useSubstream && ~useDefaultStream && ~streamsOnPool).
%   This occasion will arise if the RNGscheme was created for serial
%   computation, and given an explicit scalar RandStream object.
%   If either of these two tests fails, the function returns false;
%   however, the RNGscheme is not itself marked invalid, since it may
%   be valid with the opposite value of useParallel.
%
%   It is unnecessary to execute this function if it is known that the
%   context in which the RNGscheme was created has not changed when
%   it is now about to be used.  This may often be the case when the 
%   RNGscheme is created and used in the same toolbox command, and all
%   the code is property of a knowledgleable developer.
%   However, the validity check is wise in any kind of general purpose 
%   application or utility, and it will certainly be wise in any circumstance
%   where an RNGscheme structure is used across multiple toolbox commands.

%   Copyright 2010-2014 The MathWorks, Inc.

valid = true;

if isempty(RNGscheme)
    % Parallel or serial computation w/o RNG specification is ruled valid.
    return
end

if ~RNGscheme.valid
    valid = false;
    error(message('stats:internal:iscompatibleRNGscheme:invalidRNGscheme', RNGscheme.uuid));
end

useSubstreams    = RNGscheme.useSubstreams;
streams          = RNGscheme.streams;
useDefaultStream = RNGscheme.useDefaultStream;
streamsOnPool    = RNGscheme.streamsOnPool;

if streamsOnPool && ~useParallel
    valid = false;
    error(message('stats:internal:iscompatibleRNGscheme:NeedsParallel', RNGscheme.uuid));
    return
end

if useParallel && ~useDefaultStream && ~useSubstreams 
    % We don't create an RNGscheme with multiple streams unless
    % a parallel pool is open at the time of creation.  Therefore, the creator
    % of the RNGscheme is licensed for PCT and unguarded getParallelPoolSize
    % should be valid.  No guarantees are made for an RNGscheme struct that
    % is saved and loaded back by a user without a PCT license.
    if length(streams)>1 && length(streams) ~= internal.stats.parallel.getParallelPoolSize
        % This RNGscheme was created for a different parallel pool than the
        % one that is open now.
        RNGscheme.valid = false;
        valid = false;
        error(message('stats:internal:iscompatibleRNGscheme:PoolChanged', RNGscheme.uuid));
        return;
    end
end

end %-iscompatibleRNGscheme
