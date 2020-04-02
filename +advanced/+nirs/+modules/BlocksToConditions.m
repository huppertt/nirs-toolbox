classdef BlocksToConditions < nirs.modules.AbstractModule
%%  BlocksToConditions - Separate blocks into their own conditions
%   Eg., Block #2 of Condition #3 becomes a condition named "cond-3 ◄ block-0002"

    properties
    end
    
    methods

        function obj = BlocksToConditions( prevJob )
           obj.name = 'Separate blocks into their own conditions';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            

            for sub = 1:length(data)
                
                stimdata = data(sub).stimulus;
                conds = stimdata.keys;
                data(sub).stimulus = Dictionary;
                
                for c = 1:length(conds)
                    
                    stim = stimdata(conds{c});
                    onset = stim.onset;
                    dur = stim.dur;
                    amp = stim.amp;
                    nblock = length(onset);
                    
                    for b = 1:nblock
                        name = [conds{c} ' ◄ ' sprintf('block-%06i',b)];
                        data(sub).stimulus(name) = nirs.design.StimulusEvents(name,onset(b),dur(b),amp(b));
                    end
                    
                end
            end
        end
    end
end
