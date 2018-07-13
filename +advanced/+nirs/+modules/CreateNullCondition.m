classdef CreateNullCondition < nirs.modules.AbstractModule
    %% CreateNullCondition - Creates a condition filling out all empty time
    %
    % Options:
    %     ConditionName - Name to use for null condition (Default: 'Null')
    %
    properties
        ConditionName = 'Null';
    end
    
    methods
        function obj = CreateNullCondition( prevJob )
            if nargin > 0; obj.prevJob = prevJob; end
            obj.name = 'Craetes a condition filling out all empty time';
        end
        
        function data = runThis( obj, data )
            
            for i = 1:numel(data)
                
                X = data(i).getStimMatrix;
                
                null_time = ~any(X,2);
                delta_null = diff([0; null_time; 0]);
                onsets = find(delta_null==1);
                offsets = find(delta_null==-1);
                
                assert(length(onsets)==length(offsets),'Condition timing mismatch');
                
                if isempty(onsets)
                    continue;
                end
                
                t = data(i).time;
                Fs = data(i).Fs;
                t(end+1) = t(end)+1/Fs;
                
                stim = nirs.design.StimulusEvents(obj.ConditionName);
                stim.onset = t(onsets);
                stim.dur = t(offsets) - t(onsets);
                stim.amp = ones(size(onsets));
                
                assert(isequal(stim.getStimVector(data(i).time),null_time),'Reconstructed null vector doesn''t match original');

                data(i).stimulus(obj.ConditionName) = stim;
                
            end
        end
    end 
end