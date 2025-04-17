classdef DiscardTypes < nirs.modules.AbstractModule
    %% DiscardTypes - Removes data types (HbO/HbR/HbT)
    %
    % Options:
    %     Types - String of type to keep or cell array of types

    properties
        types = []; % char or cell array of type to keep
    end

    methods
        function obj = DiscardTypes( prevJob )
            if nargin > 0; obj.prevJob = prevJob; end
            obj.name = 'Keep Just These Types';
        end

        function data = runThis( obj, data )

            if isempty(obj.types), warning('No types specified. Skipping.'); return; end
            if ischar(obj.types), obj.types = {obj.types}; end

            for i = 1:numel(data)
                if(isa(data(i),'nirs.core.sFCStats'))
                    tbl=data(i).table;
                    lst=find(ismember(tbl.TypeDest,obj.types) | ismember(tbl.TypeOrigin,obj.types));
                    lstc=find(ismember(tbl.condition,tbl.condition{1}));
                    lsta=find(ismember(tbl(lstc,:).TypeDest,obj.types) | ismember(tbl(lstc,:).TypeOrigin,obj.types));

                    if(any(ismember(obj.types,'same')))
                        if(iscell(tbl.TypeDest))
                            lst=[lst; find(strcmp(tbl.TypeOrigin,tbl.TypeDest))];

                            lsta=[lsta; find(strcmp(tbl(lstc,:).TypeOrigin,tbl(lstc,:).TypeDest))];

                        else
                            lst=[lst; find(tbl.TypeOrigin==tbl.TypeDest)];
                            lsta=[lsta; find(tbl(lstc,:).TypeOrigin==tbl(lstc,:).TypeDest)];

                        end

                    end

                    if(any(ismember(obj.types,'diff')))
                        if(iscell(tbl.TypeDest))
                            lst=[lst; find(~strcmp(tbl.TypeOrigin,tbl.TypeDest))];


                            lsta=[lsta; find(~strcmp(tbl(lstc,:).TypeOrigin,tbl(lstc,:).TypeDest))];
                        else
                            lst=[lst; find(tbl.TypeOrigin~=tbl.TypeDest)];

                            lsta=[lsta; find(tbl(lstc,:).TypeOrigin~=tbl(lstc,:).TypeDest)];

                        end

                    end
                    data(i).probe.connections(lst,:)=[];
                    data(i).R(lst)=[];
                    if(~isempty(data(i).ZstdErr))
                        data(i).ZstdErr(lsta,:,:)=[];
                    end
               
                else
                    switch class(data)
                        case {'nirs.core.Data'}
                            tdata = data(i).probe.link.type;
                        case {'nirs.core.ChannelStats','nirs.core.ChannelFStats'}
                            tdata = data(i).variables.type;
                    end
                    tdata_probe = data(i).probe.link.type;

                    keep_type = true(size(tdata));
                    keep_type_probe = true(size(tdata_probe));
                    for j = 1:length(obj.types)
                        keep_type( strcmpi( tdata , obj.types{j} ) ) = false;
                        keep_type_probe( strcmpi( tdata_probe , obj.types{j} ) ) = false;
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

                    end
                end
            end
        end
    end
end
