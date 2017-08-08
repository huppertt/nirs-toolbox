function [useParallel,RNGscheme,poolSize] = ...
    processParallelAndStreamOptions(opt, multipleStreams)
%PROCESSPARALLELANDSTREAMOPTIONS vets and organizes parallel and stream info.
% 
%   The command line arguments together define a set of ground rules as to 
%   how random number streams are to be utilized during a sequence of computations.
%   The input argument "opt" supplies the first three of the effective
%   arguments listed below.
%
%   The effective arguments are:
%
%      useParallel        boolean      Use parfor-loops in place
%                                      of for-loops (T/F).
%      useSubstreams      boolean      Use a separate Substream for each iterate
%                                      of the loop, both in for-loops and in
%                                      parfor-loops) (T/F).
%      streams            cell array   RandStream object(s), can be empty.
%      multipleStreams    boolean      Allow multiple user-supplied streams 
%                                      to be distributed to workers in a 
%                                      parallel pool (T/F).  Some toolbox commands
%                                      will support this paradigm, others not.
%
%   The return values are:
%
%      useParallel        boolean      Use parfor loops (T/F) vs for-loops.
%      RNGscheme          struct       Contains essential and convenience 
%                                      information defining the ground rules, 
%                                      or "scheme" for RNG use.
%      poolsz             integer      The size of the parallel pool to be used
%                                      for computation.  The value is zero if
%                                      no parallel pool is open or if useParallel
%                                      is false.
%
%   Contents of the RNGscheme struct are:
%      uuid               string       Unique identifier for the RNGscheme
%      useSubstreams      boolean      As described above
%      streams            cell array   As described above
%      valid              boolean      RNGscheme is currently valid 
%      useDefaultStream   boolean      RNGscheme uses the default RandStream
%      streamsOnPool      boolean      RNGscheme deploys multiple streams on the 
%                                      parallel pool
%
%   Notes:
%      1. The purpose of the field "uuid" is to support persistence of the
%         RNGscheme and allow multiple RNGschemes to be available at the
%         same time.
%      2. The value for "uuid" is derived using the command "tempname".
%         The documentation for "tempname" indicates that, when running
%         MATLAB without a JVM, the value for "uuid" is almost
%         certain to be unique, but the guarantee is not absolute.
%         We only require uniqueness for the lifetime of the MATLAB executable
%         in which the RNGscheme is created (this is the client executable
%         when there is a parallel pool open).
%      3. The field "valid" is intended to be writeable.  All other fields
%         should be immutable, once the RNGscheme struct has been initialized.
%      4. An RNGscheme will become invalid under the following conditions:
%         a. The RNGscheme was created for parallel computation with multiple 
%            streams and an open parallel pool, and that parallel pool has since 
%            been closed.
%         b. The RNGscheme was created for parallel computation with multiple 
%            streams and an open parallel pool.  Subsequently, serial computation
%            is requested, using the RNGscheme.
%         c. The RNGscheme was created without a parallel pool, but with a
%            RandStream object supplied in the command line argument "streams".
%            Subsequently, a parallel pool is open, and parallel computation
%            is attempted, using the RNGscheme.

%   Copyright 2010-2014 The MathWorks, Inc.

if nargin<2
    multipleStreams = false;
end

[useParallel, useSubstreams, streams] = ...
    internal.stats.parallel.extractParallelAndStreamFields(opt);

% Check for valid Options parameters
[streamsOnPool,poolSize,useParallel] = parforValidateStreamOptions(useParallel, ...
                                                     useSubstreams, ...
                                                     streams, ...
                                                     multipleStreams);

% Create and initialize the return argument
phonydir = ['phony' filesep];
uuid = regexprep(tempname('phony'),phonydir,'');
RNGscheme = struct('uuid',uuid, ...
                   'useSubstreams', useSubstreams, ...
                   'streams', [], ...
                   'valid', true, ...
                   'useDefaultStream', isempty(streams), ...
                   'streamsOnPool', streamsOnPool);
RNGscheme.streams = streams;

% If using multiple streams on a parallel pool, distribute one stream to
% each of the workers, retrievable via the key "uuid".
if streamsOnPool
    % A parallel pool is open.
    % The command line supplied multiple streams.
    % We need to fan the streams out to the parallel pool.
    internal.stats.parallel.distributeToPool(uuid,streams);
end

end   % of processParforOptions()

function [streamsOnPool,poolSize,useParallel] = parforValidateStreamOptions( ...
    useParallel, useSubstreams, streams, multipleStreams ) 
%
%   This is a utility function used to support statistics functions that may
%   employ parfor-loops. The function checks that options affecting random 
%   number usage are valid in the current execution environment.
%   Invalid combinations of options cause an exception to be thrown.
%   The options to be validated are:
%
%      useParallel      boolean       Use parfor-loops in place
%                                     of for-loops (T/F).
%      useSubstreams    boolean       Use a separate Substream for each iterate
%                                     of the loop, both in for-loops and in
%                                     parfor-loops) (T/F).
%      streams          cell array    RandStream object(s), can be empty.
%      multipleStreams  boolean       Allow multiple user-supplied streams 
%                                     to be distributed to workers in a 
%                                     parallel pool (T/F).  Some toolbox commands
%                                     will support this paradigm, others not.
%
%   The return values are:
%
%      streamsOnPool    boolean       Multiple user-supplied streams are to be
%                                     used by workers in the parallel pool (T/F).
%      poolsz           integer       The size of the parallel pool to be used
%                                     for computation.  The value is zero if
%                                     no parallel pool is open or if useParallel
%                                     is false.
%
%   Factors determining validity are:
%
%      (1) Availability of the Parallel Computing Toolbox (PCT)
%          (ie, licensing and installation)
%      (2) Presence or absence of an open parallel pool
%      (3) Statistics and Machine Learning Toolbox guidelines for random
%          number stream usage in serial and parallel contexts.
%
%   The rules determining validity are:
%
%      (1) useParallel (ie, parfor) is valid w/wo PCT
%      (2) useParallel (ie, parfor) is valid w/wo a parallel pool
%          (parfor loops run in serial on the client if no parallel pool)
%      (3) streams must be empty or scalar unless all of the following hold
%             useParallel     is true
%             a parallel pool is open
%             multipleStreams is true
%             useSubstreams   is false
%          in which case
%             length(streams) must equal the parallel pool size
%          and
%             we return streamsOnPool = true
%             (otherwise, streamsOnPool = false)
%      (4) If useSubstreams is true the stream type must support Substreams.
%          The calling convention of this function is that if useSubstreams
%          is true, then the streams argument cannot be empty; a scalar
%          value must be supplied.

% devolve to serial if no parallel environment
usePool = false;
poolSize = 0;
streamsOnPool = false;

if useParallel
%     try
%         % gcp() will error if the Parallel Computing Toolbox is not
%         % installed, or if it is unable to check out a license
%         currPool = gcp();
%         if isempty(currPool)
%             poolSize = 0;
%             useParallel = false;
%         else
%             poolSize = currPool.NumWorkers;
%             usePool = true;
%         end        
%     catch
%         % If we are here, gcp() errored.  Just keep going.
%         % We have already assumed that there is no pool.
%         useParallel = false;
%     end
    poolSize = internal.stats.parallel.getParallelPoolSize('Guarded',true);
    if poolSize == 0
        useParallel = false;
        usePool = false;
    else
        usePool = true;
    end
end

if useSubstreams
    % Can only use a single stream, regardless of serial/parallel, 
    % number of workers.
    if length(streams)>1
        error(message('stats:internal:processParallelAndStreamOptions:MultipleStreamsSubstreams'));
    end
    % Make sure that this RandStream type supports Substreams.
    % Do so by seeing if it is possible to change the Substream property
    % on a stream of the same type.
    s = streams{1};
    sisterStream = RandStream(s.Type);
    try
        sisterStream.Substream = sisterStream.Substream+1;
    catch ME
        clear sisterStream;
        throw(ME);
    end
    clear sisterStream
    return
end

if ~isempty(streams)
    % A 'Streams' parameter was supplied.
    if length(streams)>1
        if ~multipleStreams
            error(message('stats:internal:processParallelAndStreamOptions:MultipleStreamsFunction'));
        elseif ~usePool
            error(message('stats:internal:processParallelAndStreamOptions:MultipleStreamsSerial'));
        end
    end
    
    if usePool && multipleStreams
        % There is an open parallel pool and the function supports
        % distributing multiple streams to the pool.
        % Enforce the condition that the number of
        % supplied RandStream objects must match the parallel pool size.
        if length(streams) ~= poolSize
            error(message('stats:internal:processParallelAndStreamOptions:MultipleStreamsPoolSize'));
        end
        % This variable records that Streams should be distributed to
        % the parallel pool, one stream for each worker.
        streamsOnPool = true;
    end
end

end % of parforValidateStreamOptions()
