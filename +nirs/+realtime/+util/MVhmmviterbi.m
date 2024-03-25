function [currentState, logP] = MVhmmviterbi(seq,tr,e,varargin)
%HMMVITERBI calculates the most probable state path for a sequence.
%   STATES = HMMVITERBI(SEQ,TRANSITIONS,EMISSIONS) given a sequence, SEQ,
%   calculates the most likely path through the Hidden Markov Model
%   specified by transition probability matrix, TRANSITIONS, and emission
%   probability matrix, EMISSIONS. TRANSITIONS(I,J) is the probability of
%   transition from state I to state J. EMISSIONS(K,L) is the probability
%   that symbol L is emitted from state K.
%
%   HMMVITERBI(...,'SYMBOLS',SYMBOLS) allows you to specify the symbols
%   that are emitted. SYMBOLS can be a numeric array or a string/cell array
%   of the names of the symbols.  The default symbols are integers 1
%   through N, where N is the number of possible emissions.
%
%   HMMVITERBI(...,'STATENAMES',STATENAMES) allows you to specify the names
%   of the states. STATENAMES can be a numeric array or a cell array of the
%   names of the states. The default statenames are 1 through M, where M is
%   the number of states.
%
%   This function always starts the model in state 1 and then makes a
%   transition to the first step using the probabilities in the first row
%   of the transition matrix. So in the example given below, the first
%   element of the output states will be 1 with probability 0.95 and 2 with
%   probability .05.
%
%   Examples:
%
%       tr = [0.95,0.05; 0.10,0.90];
%
%       e = [1/6,  1/6,  1/6,  1/6,  1/6,  1/6; ...
%            1/10, 1/10, 1/10, 1/10, 1/10, 1/2;];
%
%       [seq, states] = hmmgenerate(100,tr,e); 
%       estimatedStates = hmmviterbi(seq,tr,e);
%
%       [seq, states] = ...
%       hmmgenerate(100,tr,e,'Statenames',{'fair';'loaded'});
%       estimatesStates = ...
%       hmmviterbi(seq,tr,e,'Statenames',{'fair';'loaded'});
%
%   See also HMMGENERATE, HMMDECODE, HMMESTIMATE, HMMTRAIN.

%   Reference: Biological Sequence Analysis, Durbin, Eddy, Krogh, and
%   Mitchison, Cambridge University Press, 1998.

%   Copyright 1993-2011 The MathWorks, Inc.


% tr must be square

if nargin > 0
    if isstring(seq)
        seq = cellstr(seq);
    end
end

if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

numStates = size(tr,1);
checkTr = size(tr,2);
if checkTr ~= numStates
    error(message('stats:hmmviterbi:BadTransitions'));
end

% number of rows of e must be same as number of states

checkE = length(e);
if checkE ~= numStates
    error(message('stats:hmmviterbi:InputSizeMismatch'));
end

customStatenames = false;

% deal with options
if nargin > 3
    okargs = {'symbols','statenames'};
    [symbols,statenames] = ...
        internal.stats.parseArgs(okargs, {[] []}, varargin{:});
    
    if ~isempty(symbols)
        numSymbolNames = numel(symbols);
        if ~isvector(symbols) || numSymbolNames ~= numEmissions
            error(message('stats:hmmviterbi:BadSymbols'));
        end
        [~, seq]  = ismember(seq,symbols);
        if any(seq(:)==0)
            error(message('stats:hmmviterbi:MissingSymbol'));
        end
    end
    if ~isempty(statenames)
        numStateNames = length(statenames);
        if numStateNames ~= numStates
            error(message('stats:hmmviterbi:BadStateNames'));
        end
        customStatenames = true;
    end
end


% work in log space to avoid numerical issues
L = size(seq,2);

currentState = zeros(1,L);
if L == 0
    return
end
logTR = log(tr);


% allocate space
pTR = zeros(numStates,L);
% assumption is that model is in state 1 at step 0
v = -Inf(numStates,1);
v(1,1) = 0;
vOld = v;
Vpath= zeros(numStates,L);
% loop through the model
for count = 1:L
    for state = 1:numStates
        % for each state we calculate
        % v(state) = e(state,seq(count))* max_k(vOld(:)*tr(k,state));
        bestVal = -inf;
        bestPTR = 0;
        % use a loop to avoid lots of calls to max
        for inner = 1:numStates 
            val = vOld(inner) + logTR(inner,state);
            if val > bestVal
                bestVal = val;
                bestPTR = inner;
            end
        end
        % save the best transition information for later backtracking
        pTR(state,count) = bestPTR;
        % update v
        logE = log(e{state}.pdf(seq(:,count)'));
       % -1/2*seq(:,count)'*inv(squeeze(e(state,:,:)))*seq(:,count) + log(1/(sqrt(2*pi)*sqrt(norm(squeeze(e(state,:,:))))));
        v(state) = logE + bestVal;
    end
    vOld = v;
    Vpath(:,count)=v;
end

% decide which of the final states is post probable
[logP, finalState] = max(v);

% Now back trace through the model
currentState(L) = finalState;
for count = L-1:-1:1
    currentState(count) = pTR(currentState(count+1),count+1);
    if currentState(count) == 0
        error(message('stats:hmmviterbi:ZeroTransitionProbability', currentState( count + 1 )));
    end
end
if customStatenames
    currentState = reshape(statenames(currentState),1,L);
end


