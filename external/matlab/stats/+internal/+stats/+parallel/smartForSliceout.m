function varargout = smartForSliceout( niter, loopbody, useParallel, RNGscheme)
%SMARTFORSLICEOUT Adaptive dispatch of for-loops or parfor-loops.
%
%   SMARTFORSLICEOUT provides dispatching of for-loops or parfor-loops on
%   behalf of Statistics and Machine Learning Toolbox commands that support
%   both serial and parallel modes of computation.
%
%   VARARGOUT = SMARTFORSLICEOUT(NITER,LOOPBODY,USEPARALLEL) executes a FOR-loop 
%   if USEPARALLEL is FALSE, or a PARFOR-loop if USEPARALLEL is TRUE.
%   In either case, the index range for the loop is 1:NITER,
%   and each iteration of the loop performs one invocation of
%   the function handle LOOPBODY.  See below for details about LOOPBODY.
%   NITER is a positive integer.
%   LOOPBODY supports a calling signature as described below.
%   USEPARALLEL is logical TRUE or FALSE.
%
%   VARARGOUT = SMARTFORSLICEOUT(...,RNGSCHEME) specifies a struct RNGSCHEME, which
%   contains information about how the LOOPBODY function is to generate 
%   random numbers.  Normally RNGSCHEME is obtained by calling the utility
%   function INTERNAL.STATS.PARALLEL.PROCESSPARALLELANDSTREAMOPTIONS(OPT),
%   where OPT is a user-supplied command line argument to the toolbox
%   command that calls SMARTFORSLICEOUT.  RNGSCHEME may be absent or empty, 
%   if LOOPBODY will not use random number generators, or will only 
%   use the MATLAB default random number stream(s).
%
%   The function LOOPBODY must support the calling signature
%   VARARGOUT = LOOPBODY(ITER,S), where ITER is the index of the
%   for-loop or parfor-loop executing within SMARTFORSLICEOUT, 
%   and S is either empty or a RandStream object.
%   If S is empty, it tells LOOPBODY that no user-supplied RANDSTREAM
%   object was specified, in which case, if LOOPBODY needs to generate
%   random numbers, it should use the MATLAB default random stream.
%
%   The contents of the struct RNGSCHEME are:
%   
%       uuid               string       Unique identifier for the RNGscheme 
%       useSubstreams      boolean      If 'true', use a separate Substream 
%                                       for each iterate of the loop
%       streams            cell array   The RandStream object(s) to be used 
%                                       for random number generation during the
%                                       invocation of SMARTFORSLICEOUT. Can be empty. 
%       valid              boolean      RNGSCHEME is currently valid (T/F)
%       useDefaultStream   boolean      RNGSCHEME uses default RANDSTREAM(s) (T/F)
%       streamsOnPool      boolean      RNGscheme deploys multiple streams 
%                                       on a parallel pool (T/F)
%
%   SMARTFORSLICEOUT is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.

%   Copyright 2010-2014 The MathWorks, Inc.

%-------- Initial sanity check --------

% For efficiency we skip these, but can provide checks if warranted:
%     1. niter is a positive integer
%     2. loopbody is a function handle
%     3. useParallel is logical

%-------- If present, check and pre-process 'RNGscheme' --------

if nargin < 4 
    RNGscheme = [];
end

%-------- If present, check and pre-process 'RNGscheme' --------

if ~isempty(RNGscheme)
    % Abort if we know that the RNG inputs will not work with the
    % requested computation mode.
    internal.stats.parallel.iscompatibleRNGscheme(useParallel,RNGscheme);
end

[streamsOnPool, useSubstreams, S, uuid] = ...
    internal.stats.parallel.unpackRNGscheme(RNGscheme);

if useSubstreams
    initialSubstream = internal.stats.parallel.freshSubstream(S);
else
    % This value will not be used but the definition is 
    % necessary for parfor static analysis.
    initialSubstream = 0;
end

usePool = useParallel && streamsOnPool;
    
%-------- Additional preamble --------

% This is a limit based on the need to declare sliced output
% variable with concrete names, when operating within a parfor
% loop.  If necessary for future use, increase this limit, and
% add more assignments to those for o, o1, o2, etc,
% in the assignment to sliced outputs when useParallel is true.
maxSlicedOutputs = 6;

% Confirm that we are able to handle the number of sliced output arguments
nSlicedOutputs = nargout;
if nSlicedOutputs > maxSlicedOutputs 
    error(message('stats:smartFor:Internal'));
end

% Initialize default shape parameters for reduction operation
% within parfor.
shapeParameters = cell(1,nSlicedOutputs);

% This definition keeps mlint happy.
collectShape = @internal.stats.parallel.collectShape;

%-------- Now embark on the loops --------

try

    if useParallel

        parfor iter = 1:niter
            
            %---- Get the RandStream to use for this iterate, if any ----
            T = internal.stats.parallel.prepareStream( ...
                    iter,initialSubstream,S,usePool,useSubstreams,uuid);

            %---- Run the loop body ---
            %
            % We add a fixed lead dimension so that parfor allows us
            % to return the output of loopbody into a linearly indexed
            % cell array.
            loopbodyout = cell(1,2);
            loopbodyout{1} = cell(1,nSlicedOutputs);
            [loopbodyout{1}{:}] = feval(loopbody, iter, T);
            
            slice = loopbodyout{1};

            %---- Handle sliced outputs ----
            %
            % Collect the shape parameters of the outputs
            if iter == 1
                shapeParameters = collectShape(shapeParameters,slice);
            end
            
            % Pad slice{} to length maxSlicedOutputs
            if nSlicedOutputs < maxSlicedOutputs
                slice(nSlicedOutputs+1:maxSlicedOutputs) = {0};
            end
            
            
            % Assign to sliced output variables.
            % The indexing is always 2-D and sliced by column,
            % but if the entries in slice are scalar, then the
            % assignment to o* is scalar, and the resultant
            % sliced output variable is 1-D, 1 x niter.
            o(:,iter) = slice{1}(:);
            o1(:,iter) = slice{2}(:);
            o2(:,iter) = slice{3}(:);
            o3(:,iter) = slice{4}(:);
            o4(:,iter) = slice{5}(:);
            o5(:,iter) = slice{6}(:);
            % Add assignments here for o6, o7, etc, if it is
            % necessary to increase maxSlicedOutputs.
            % ...
            
        end %-of parfor loop

    else

        %-------- SERIAL COMPUTATION --------

        % This will hold the return values from a single loop iteration.
        slice = cell(1,nSlicedOutputs);
        
        for iter = 1:niter
            if useSubstreams
                S.Substream = iter + initialSubstream - 1;
            end

            [slice{:}] = loopbody(iter, S);

            if iter==1
                for i=1:nSlicedOutputs
                    % Preallocate storage depending on the sizing and
                    % type of the first set of sliced outputs.
                    % Assumption: every accumulated variable is "uniform",
                    % ie, each column has the same length.
                    varargout{i}(:,niter) = slice{i}(:);
                    
                    % Retain the shape parameters of all the sliced
                    % outputs of loopbody.
                    shapeParameters{i} = size(slice{i});
                end                    
            end
            
            % Assign the computed slices to the real output variables
            for i=1:nSlicedOutputs
                varargout{i}(:,iter) = slice{i}(:);
            end
        end

    end
    
    %-------- Post-loop coordination --------
    if ~isempty(RNGscheme)
        internal.stats.parallel.reconcileStreamsAfterLoop( ...
            RNGscheme, initialSubstream, niter);
    end
    
    if useParallel
        varargout = {o o1 o2 o3 o4 o5};
        varargout = varargout(1:nSlicedOutputs);
    end

    % Reshape the sliced outputs
    for i=1:nSlicedOutputs
        if shapeParameters{i}(1) == 1 && length(shapeParameters{i}) == 2
            % We got row vectors or scalars and we transpose
            % to return the result shaped as [niter m].
            varargout{i} = varargout{i}';
        elseif length(shapeParameters{i}) == 2 && shapeParameters{i}(2) == 1
            % Loopbody returned a column vector dimensioned [m 1]. 
            % We have already accumulated iterates the way
            % we want them (ie, result is shaped [m niter]).
            % Do nothing.
        else
            % We stack along the rightmost dimension.
            % If the result of one iteration was shaped [m n p]
            % we will return an output of shape [m n p niter].
            % This will be true even if m or n or p is one.
            varargout{i} = reshape(varargout{i}, [shapeParameters{i} niter]);
        end
    end

catch ME
    % If Substreams were used, reset the Substream as if we had
    % run all niter iterates.
    % Note that as this try/catch is currently positioned,
    % we will reset the Substream if an error occurs in the
    % parfor/for loop, but not if an error occurs in the 
    % preamble part of the function, before any iterates are
    % performed (in which case the Substream was not altered
    % within this routine).  We could consider resetting the 
    % Substream in the same way regardless of when/where an error
    % occurs.  However, since errors can potentially occur
    % in processing the RNGscheme argument, we may not always 
    % have a valid definition for the random number stream (S).
    if useSubstreams
        S.Substream = initialSubstream + niter;
    end
    rethrow(ME);
end

end %-smartFor

