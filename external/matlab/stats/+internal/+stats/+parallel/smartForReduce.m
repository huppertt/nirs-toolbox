function reductionValue = smartForReduce( niter, loopbody, useParallel, ...
                                    RNGscheme, reductionVariables)
%SMARTFORREDUCE Adaptive dispatch of for-loops or parfor-loops.
%
%   SMARTFORREDUCE provides dispatching of for-loops or parfor-loops on
%   behalf of Statistics and Machine Learning Toolbox commands that support
%   both serial and parallel modes of computation.
%
%   VARARGOUT = ...
%       SMARTFORREDUCE(NITER,LOOPBODY,USEPARALLEL,RNGSCHEME,REDUCTIONVARIABLES)
%   
%   REDUCTIONVARIABLES identifies reduction variables and specifies how 
%   they are to be calculated.  REDUCTIONVARIABLES may have one of two forms.  
%   The canonical form is a 2-element cell array.  In this form, the first 
%   element of the cell array is a function handle, REDUCTFUN.
%   REDUCTFUN is of the form VAL = REDUCTFUN(VAL,UPDATE), 
%   and is commutative and associative. The second element of the cell
%   array is an initial value for the reduction variable, VAL, so that
%   the first invocation of REDUCTFUN will have initialized operands.  
%   The second form for REDUCTIONVARIABLES is a keyword string, which can 
%   take the values 'argmin' or 'argmax'.  If the keyword form is specified, 
%   SMARTFORREDUCE will supply the appropriate REDUCTFUN and initial value.  
%   If REDUCTIONVARIABLES is present, and not empty, the function LOOPBODY 
%   must supply the reduction variable value calculated on each invocation
%   as the last, or rightmost, of its return arguments.  If the keyword
%   'argmin' or 'argmax' is specified, the reduction variable value
%   must be a vector cell array of length 2 or greater.  The first two 
%   elements of the cell array must contain a numeric value that can be 
%   supplied as an operand to relational operators (lt or gt, respectively).
%   The second element is used as a tie-breaker in case of numerical 
%   equality between reduction variable values in two different loop
%   iterations.  The recommend practice is that LOOPBODY place the loop 
%   iteration in the second index of the cell array that represents the 
%   reduction variable value.  This practice will cause SMARTFORREDUCE
%   to return the result of the least iteration (with 'argmin') or 
%   greatest (with 'argmax'), in case of ties.  This result, in turn,
%   permits SMARTFORREDUCE to return the same result in both serial
%   and parallel computation.
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
%                                       invocation of SMARTFORREDUCE. Can be empty. 
%       valid              boolean      RNGSCHEME is currently valid (T/F)
%       useDefaultStream   boolean      RNGSCHEME uses default RANDSTREAM(s) (T/F)
%       streamsOnPool      boolean      RNGscheme deploys multiple streams 
%                                       on a parallel pool (T/F)
%
%   SMARTFORREDUCE is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.

%   Copyright 2010-2014 The MathWorks, Inc.

    %-------- Initial sanity check --------
    
    % For efficiency we skip these, but can provide checks if warranted:
    %     1. niter is a positive integer
    %     2. loopbody is a function handle
    %     3. useParallel is logical

if nargin < 5 || isempty(reductionVariables)
    % This function should not be called unless reduction variables
    % are to be calculated.
    error(message('stats:parallel:smartForReduce'));
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

%-------- Pre-process 'reductionVariables' --------
    
[reductionOperator,reductionValue] = ...
     internal.stats.parallel.processReductionVariableArgument(reductionVariables);

%-------- Now embark on the loops --------

try

    if useParallel
        
        parfor iter = 1:niter
            
            %---- Get the RandStream to use for this iterate, if any ----
            T = internal.stats.parallel.prepareStream( ...
                iter,initialSubstream,S,usePool,useSubstreams,uuid);

            %---- Run the loop body ---
            reduce = feval(loopbody, iter, T);
            
            %---- Update any reduction variables ----
            reductionValue = reductionOperator(reductionValue, reduce);
            
        end %-of parfor loop

    else

        %-------- SERIAL COMPUTATION --------

        for iter = 1:niter

            if useSubstreams
                S.Substream = iter + initialSubstream - 1;
            end
            
            reduce = loopbody(iter, S);
            
            reductionValue = reductionOperator(reductionValue, reduce);
        end
        
    end
    
    %-------- Post-loop coordination --------
    if ~isempty(RNGscheme)
        internal.stats.parallel.reconcileStreamsAfterLoop( ...
            RNGscheme, initialSubstream, niter);
    end

catch ME
    % If Substreams were used, reset the Substream as if we had
    % run all niter iterates.
    % Note that as this try/catch is currently positioned,
    % we will reset the Substream if an error occurrs in the
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

