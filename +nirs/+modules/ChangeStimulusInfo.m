classdef ChangeStimulusInfo < nirs.modules.AbstractModule
%%  ChangeStimulusInfo - Changes sitmulus info to data given a table.
% 
% Options: 
%     ChangeTable - table of the stim info to replace in the data (use NaN
%     to keep old values; e.g. to chnage only the duration and not the
%     timing)

    properties
        ChangeTable =table();   % a matlab table with the stim info
    end
    
    methods

        function obj = ChangeStimulusInfo( prevJob )
           obj.name = 'Change Stimulus Information';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            % columns of the table that arent varToMatch
           oldstiminfo = nirs.createStimulusTable(data);
           newstiminfo = obj.ChangeTable;
           
           if(height(newstiminfo)==1 & isnan(newstiminfo.FileIdx))
               % assume you ment to use this for all files 
               % replicate the table
               newstiminfo=repmat(newstiminfo,length(data),1);
               newstiminfo.FileIdx=[1:length(data)]';
           end
           
           VarNames=newstiminfo.Properties.VariableNames;
           VarNames={VarNames{~ismember(VarNames,'FileIdx')}};
           for idx=1:height(newstiminfo)
               fileIdx=newstiminfo.FileIdx(idx);
               
                for idx2=1:length(VarNames)
                    stim = data(fileIdx).stimulus(VarNames{idx2});
                    newstim=newstiminfo.(VarNames{idx2})(idx);
                    if(~isa(stim,class(newstim)))
                        stim=newstim;
                    else
                        if(~isnan(newstim.onset))
                            stim.onset=newstim.onset;
                        end
                        if(~isnan(newstim.dur))
                            if(length(newstim.dur)==1)
                                newstim.dur=newstim.dur*ones(size(stim.onset));
                            end
                            stim.dur=newstim.dur;
                        end
                        if(~isnan(newstim.amp))
                            if(length(newstim.amp)==1)
                                newstim.amp=newstim.amp*ones(size(stim.onset));
                            end
                            stim.amp=newstim.amp;
                        end
                    end
                    
                    data(fileIdx).stimulus(VarNames{idx2})=stim;    
                end
           end
            
            
        end
    end
    
end