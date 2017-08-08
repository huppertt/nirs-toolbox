classdef KeepTypes < nirs.modules.AbstractModule
    %% KeepTypes - Removes all data types (HbO/HbR/HbT) except those specified.
    %
    % Options:
    %     Types - String of type to keep or cell array of types

    properties
        types = []; % char or cell array of type to keep
    end
    
    methods
        function obj = KeepTypes( prevJob )
            if nargin > 0; obj.prevJob = prevJob; end
            obj.name = 'Keep Just These Types';
        end
        
        function data = runThis( obj, data )
            
            if isempty(obj.types), warning('No types specified. Skipping.'); return; end
            if ischar(obj.types), obj.types = {obj.types}; end
            
            for i = 1:length(data)
                
                switch class(data)
                    case {'nirs.core.Data','nirs.core.sFCStats'}
                        tdata = data(i).probe.link.type;
                    case {'nirs.core.ChannelStats','nirs.core.ChannelFStats'}
                        tdata = data(i).variables.type;
                end
                tdata_probe = data(i).probe.link.type;
                
                keep_type = false(size(tdata));
                keep_type_probe = false(size(tdata_probe));
                for j = 1:length(obj.types)
                    keep_type( strcmpi( tdata , obj.types{j} ) ) = true;
                    keep_type_probe( strcmpi( tdata_probe , obj.types{j} ) ) = true;
                end

                data(i).probe.link = data(i).probe.link(keep_type_probe,:);
                
                switch class(data)
                    case {'nirs.core.Data'}
                        data(i).data = data(i).data(:,keep_type);
                    case {'nirs.core.ChannelStats'}
                        data(i).variables = data(i).variables(keep_type,:);
                        data(i).beta = data(i).beta(keep_type);
                        data(i).covb = data(i).covb(keep_type,keep_type);
                    case {'nirs.core.ChannelFStats'}
                        data(i).variables = data(i).variables(keep_type,:);
                        data(i).F = data(i).F(keep_type);
                    case {'nirs.core.sFCStats'}
                        data(i).R = data(i).R(keep_type,keep_type,:);
                        if(~isempty(data(i).ZstdErr))
                            data(i).ZstdErr = data(i).ZstdErr(keep_type,keep_type,:);
                        end
                end
            end
        end
    end
end
