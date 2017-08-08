function varargout = smartFor( niter, loopbody, useParallel, ...
                                     RNGscheme, reductionVariables, varargin)
%SMARTFOR Adaptive dispatch of for-loops or parfor-loops.
%
%   SMARTFOR provides dispatching of for-loops or parfor-loops on behalf of
%   Statistics and Machine Learning Toolbox commands that support both
%   serial and parallel modes of computation.
%
%   VARARGOUT = SMARTFOR(NITER,LOOPBODY,USEPARALLEL) executes a FOR-loop 
%   if USEPARALLEL is FALSE, or a PARFOR-loop if USEPARALLEL is TRUE.
%   In either case, the index range for the loop is 1:NITER,
%   and each iteration of the loop performs one invocation of
%   the function handle LOOPBODY.  NITER is a positive integer.
%   LOOPBODY supports a calling signature as described below.
%   USEPARALLEL is logical TRUE or FALSE.
%
%   VARARGOUT = SMARTFOR(...,RNGSCHEME) specifies a struct RNGSCHEME, which
%   contains information about how the LOOPBODY function is to generate 
%   random numbers.  Normally RNGSCHEME is obtained by calling the utility
%   function INTERNAL.STATS.PARALLEL.PROCESSPARALLELANDSTREAMOPTIONS(OPT),
%   where OPT is a user-supplied command line argument to the toolbox
%   command that calls SMARTFOR.  RNGSCHEME may be absent or empty, 
%   if LOOPBODY will not use random number generators, or will only 
%   use the MATLAB default random number stream(s).
%
%   VARARGOUT = SMARTFOR(...,RNGSCHEME,REDUCTIONVARIABLES) identifies 
%   reduction variables and specifies how they are to be calculated.
%   REDUCTIONVARIABLES may be absent or empty if the call to SMARTFOR 
%   does not involve any reduction variables. REDUCTIONVARIABLES may
%   have one of two forms.  The canonical form is a 2-element cell array,
%   in which the first element is a function handle REDUCTFUN of the form,
%   VAL = REDUCTFUN(VAL,UPDATE), and REDUCTFUN is commutative and 
%   associative, and the second element is an initial value for the 
%   reduction variable (eg, VAL), so that the first invocation of
%   REDUCTFUN will have initialized operands.  The second form is a 
%   keyword string, which can take the values 'argmin' or 'argmax'.  
%   If the keyword form is specified, SMARTFOR will supply the appropriate
%   REDUCTFUN and initial value.  If REDUCTIONVARIABLES is present,
%   and not empty, the function LOOPBODY must supply the reduction
%   variable value calculated on each invocation as the last, or
%   rightmost, of its return arguments.  If the keyword 'argmin' or
%   'argmax' is specified for REDUCTIONVARIABLES, the reduction variable 
%   value that LOOPBODY returns must be a cell array, and the 
%   first element of the cell array must be a numeric value that 
%   can be supplied as an operand to MIN or MAX, respectively.
%
%   VARARGOUT = SMARTFOR(...,RNGSCHEME,REDUCTIONVARIABLES,VARARGIN) 
%   specifies additional arguments to be passed by SMARTFOR to invocations
%   of LOOPBODY.  The arguments are passed in a sliced manner.  That is, 
%   for each element V in the list VARARGIN, V(:,ITER) is passed to LOOPBODY
%   on iteration number ITER of the loop dispatched by SMARTFOR.
%   If VARARGIN is present, RNGSCHEME and REDUCTIONVARIABLES must be
%   present, though they may be empty.
%   
%   The function LOOPBODY must support the calling signature
%   VARARGOUT = LOOPBODY(ITER,S,VARARGIN), where ITER is the index of the
%   for-loop or parfor-loop executing within SMARTFOR, S is either empty
%   or a RandStream object, and VARARGIN are additional input arguments 
%   that have been sliced and passed to LOOPBODY by SMARTFOR.
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
%                                       invocation of SMARTFOR. Can be empty. 
%       valid              boolean      RNGSCHEME is currently valid (T/F)
%       useDefaultStream   boolean      RNGSCHEME uses default RANDSTREAM(s) (T/F)
%       streamsOnPool      boolean      RNGscheme deploys multiple streams 
%                                       on a parallel pool (T/F)
%
%   See also SMARTFORREDUCE, SMARTFORSLICEOUT.
%
%   SMARTFOR is an internal utility and is not meant for
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

if ~isempty(RNGscheme)
% Abort if we know that the RNG inputs will not work with the
% requested computation mode.
internal.stats.parallel.iscompatibleRNGscheme(useParallel,RNGscheme)
end

[streamsOnPool, useSubstreams, S, uuid] = ...
    internal.stats.parallel.unpackRNGscheme(RNGscheme);

usePool = useParallel && streamsOnPool;

%-------- If present, pre-process 'reductionVariables' --------
    
reducedOutputs = false;
if nargin > 4 && ~isempty(reductionVariables)
    reducedOutputs = true;
    [reductionOperator,initialReductionValue] = ...
         internal.stats.parallel.processReductionVariableArgument(reductionVariables);
else
    reducedOutputs = false;
    % we need these definitions for parfor static analysis
    reductionOperator = @(a,b) a;
    initialReductionValue = 0;
end

%-------- If present, pre-process 'sliceableInputs' --------
    
if ~isempty(varargin)
    sliceableInputs = varargin;
else
    sliceableInputs = {};
end    
nSlicedInputs = length(sliceableInputs);

%-------- Additional preamble --------

% This is a limit based on the need to declare sliced output
% variable with concrete names, when operating within a parfor
% loop.  If necessary for future use, increase this limit, and
% add more assignments to those for o, o1, o2, etc,
% in the assignment to sliced outputs when useParallel is true.
maxSlicedOutputs = 6;

% Calculate the number of sliced (vs. reduced) output arguments
% and confirm that we are able to handle this many.
nSlicedOutputs = nargout - reducedOutputs;
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
    
    if useSubstreams
        initialSubstream = internal.stats.parallel.freshSubstream(S);
    else
        % This value will not be used but the definition is 
        % necessary for parfor static analysis.
        initialSubstream = 0;
    end

    if useParallel
        
        % nargout is reassigned to 1 in a parfor loop, so retain
        % the real value.
        nout = nargout;

        parfor iter = 1:niter
            
            %---- Get the RandStream to use for this iterate, if any ----
            T = internal.stats.parallel.prepareStream( ...
                    iter,initialSubstream,S,usePool,useSubstreams,uuid);

            %---- Prepare sliced inputs for @loopbody, if any ----
            if nSlicedInputs > 0
                % The extra dimension that we provide for "si"
                % lets us signal to parfor that the assignment in the
                % for-loop is not to a sliced output variable.
                si = cell(2,nSlicedInputs);
                for i=1:nSlicedInputs
                    ss = sliceableInputs{i};
                    si{1}{i} = ss(:,iter);
                end
                % And now we strip the protective blanket.
                slicedInputs = si{1};
            else
                slicedInputs = {};
            end

            %---- Run the loop body ---
            %
            % We add a fixed lead dimension so that parfor allows
            % to return the output of loopbody into a linear indexed
            % cell array.
            loopbodyout = cell(1,2);
            loopbodyout{1} = cell(1,nout);
            [loopbodyout{1}{:}] = feval(loopbody, iter, T, slicedInputs{:});
            
            %---- Separate the sliced and reduction outputs ----
            slice = loopbodyout{1}(1:nSlicedOutputs);
            if reducedOutputs
                reduce = loopbodyout(end);
            end

            %---- Handle sliced outputs, if any ----
            if nSlicedOutputs > 0
                
                % Collect the shape parameters of the outputs
                if iter == 1
                    shapeParameters = collectShape(shapeParameters, slice);
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
                
            end
            
            %---- Update any reduction variables ----
            if reducedOutputs
                initialReductionValue = ...
                    reductionOperator(initialReductionValue, reduce);
            end
            
        end %-of parfor loop

    else

        %-------- SERIAL COMPUTATION --------

        % This will hold the return values from a single loop iteration.
        loopbodyout = cell(1,nargout);
        
        if nSlicedInputs > 0
            slicedInputs = cell(1,nSlicedInputs);
        else
            slicedInputs = {};
        end
        
        for iter = 1:niter
            if useSubstreams
                S.Substream = iter + initialSubstream - 1;
            end
            
            % Extract any sliced inputs for @loopbody.
            for i=1:nSlicedInputs
                slicedInputs{i} = sliceableInputs{i}(:,iter);
            end
            
            [loopbodyout{:}] = loopbody(iter, S, slicedInputs{:});
            
            slice = loopbodyout(1:nSlicedOutputs);
            if reducedOutputs
                reduce = loopbodyout{end};
            end

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
            
            % Update any reduction variables
            if reducedOutputs
                initialReductionValue = ...
                    reductionOperator(initialReductionValue, reduce);
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
    if nSlicedOutputs > 0
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
    end

    % Assign accumulated reduction variables to varargout.
    if reducedOutputs
        varargout{nSlicedOutputs + 1} = initialReductionValue;
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
