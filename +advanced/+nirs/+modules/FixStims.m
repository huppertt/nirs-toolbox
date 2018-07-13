classdef FixStims < nirs.modules.AbstractModule
%% FixStims - Modify onset/duration/amplitude of stims
% 
% Options:
%     listOfChanges - % n x 4 cell array of changes 
%                     % (condition, onset shift, duration, amplitude)
%     
% Example:
%     j = nirs.modules.FixStims();
%     j.listOfChanges = { ...
%         'GoNoGo', -4, 30, []; ... % Subtract 4 seconds from onset, set duration to 30 seconds
%         '2-back', [], 60, 2; };   % Set duration to 60 seconds, amplitude to 2
%
%     j.run( raw );
%     
    
    properties
        listOfChanges = {}; % n x 4 cell array of changes
    end
    
    methods
        function obj = FixStims( prevJob )
           obj.name = 'Fix stimulus onsets/durations/amplitudes';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            
            if size(obj.listOfChanges,2)<4, obj.listOfChanges{1,4} = []; end
            
            assert(isa(data,'nirs.core.Data'),'Intended for nirs.core.Data objects');
            
            for i = 1:numel(data)
                
                conds = data(i).stimulus.keys;
                for j = 1:size(obj.listOfChanges,1)
                    
                    key = obj.listOfChanges{j,1};
                    ons = obj.listOfChanges{j,2};
                    dur = obj.listOfChanges{j,3};
                    amp = obj.listOfChanges{j,4};
                    
                    match_inds = find(strcmpi( conds , key ));
                    
                    for k = 1:length(match_inds)
                    
                        stim = data(i).stimulus(conds{match_inds(k)});
                        if ~isempty(ons) && ~isnan(ons)
                            stim.onset = stim.onset + ons;
                        end
                        if ~isempty(dur) && ~isnan(dur)
                            stim.dur(:) = dur;
                        end
                        if ~isempty(amp) && ~isnan(amp)
                            stim.amp(:) = amp;
                        end
                        
                        data(i).stimulus(conds{match_inds(k)}) = stim;
                    end
                
                end
            
            end
            
        end
    end
    
end

