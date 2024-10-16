classdef HiddenMarkov < nirs.realtime.modules.AbstractModule
    % real-time implementation of a band-pass filter
    
    properties
       emissions
       transitions
       type='viterbi';  % Only one coded so far
       state_labels_key; % the enum labels of the states
    end
    
    properties(Hidden=true)
        previous_state=[];;
    end    
    
    methods
        function obj=HiddenMarkov(prevJob)
            obj.name='RT-HMM';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function obj=resetThis(obj)
            obj.previousstate=[];
        end
        
        
        function [dOut,t,probe,stimulus] = updateThis(obj,d,t,probe,stimulus)

            if(isempty(obj.transitions) || isempty(obj.emissions))
                error('HMM must be trained first');
            end

            logTR=log(obj.transitions);
            numStates=size(logTR,1);

            if(isempty(obj.previous_state))
                obj.previous_state=-Inf(numStates,1);
                obj.previous_state(1,1)=0;  % Start in the first state
            end

            L=size(d,2);
            
            % allocate space
            pTR = zeros(numStates,L);
            Vpath= zeros(numStates,L);

            for count=1:L
                for state = 1:numStates
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
                    logE = log(obj.emissions{state}.pdf(d(:,count)'));
                    % -1/2*seq(:,count)'*inv(squeeze(e(state,:,:)))*seq(:,count) + log(1/(sqrt(2*pi)*sqrt(norm(squeeze(e(state,:,:))))));
                    v(state) = logE + bestVal;
                end
                obj.previous_state=v;
                Vpath(:,count)=v;
            end
            
            currentState = zeros(1,L);
            % Now back trace through the model

            % decide which of the final states is post probable
            [logP, finalState] = max(v);

            currentState(L) = finalState;
            for count = L-1:-1:1
                currentState(count) = pTR(currentState(count+1),count+1);
            end

            dout=currentState;
            % or
            %[~,dout]=max(Vpath,[],1);

            if(~isempty(obj.state_labels_key))
                dout=reshape(obj.state_labels_key(dout),1,L);
            end

        end
            
    end
    
end