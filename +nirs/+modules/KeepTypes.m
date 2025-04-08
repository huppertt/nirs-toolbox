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
            
            for i = 1:numel(data)
                if(isa(data(i),'nirs.core.sFCStats'))
                    tbl=data(i).table;
                    if(length(obj.types)==1 && strcmp(obj.types{1},'same'))
                        obj.types=unique({tbl.TypeDest{:}; tbl.TypeOrigin{:}});
                        obj.types{end+1}='same';
                    end
                    if(length(obj.types)==1 && strcmp(obj.types{1},'diff'))
                        obj.types=unique({tbl.TypeDest{:}; tbl.TypeOrigin{:}});
                        obj.types{end+1}='diff';
                    end
                    lst=find(ismember(tbl.TypeDest,obj.types) & ismember(tbl.TypeOrigin,obj.types));
                    lstc=find(ismember(tbl.condition,tbl.condition{1}));
                    lsta=find(ismember(tbl(lstc,:).TypeDest,obj.types) & ismember(tbl(lstc,:).TypeOrigin,obj.types));
                    
                    if(any(ismember(obj.types,'same')))
                        if(iscell(tbl.TypeDest))
                            lst2=find(~strcmp(tbl(lst,:).TypeOrigin,tbl(lst,:).TypeDest));
                            lst(lst2)=[];
                            
                            lst2=find(~strcmp(tbl(lsta,:).TypeOrigin,tbl(lsta,:).TypeDest));
                            lsta(lst2)=[];
                        else
                            lst2=find(tbl(lst,:).TypeOrigin~=tbl(lst,:).TypeDest);
                            lst(lst2)=[];

                            lst2=find(tbl(lstc(lsta),:).TypeOrigin~=tbl(lstc(lsta),:).TypeDest);
                            lsta(lst2)=[];

                        end

                    end

                    if(any(ismember(obj.types,'diff')))
                        if(iscell(tbl.TypeDest))
                            lst2=find(strcmp(tbl(lst,:).TypeOrigin,tbl(lst,:).TypeDest));
                            lst(lst2)=[];
                            
                            lst2=find(strcmp(tbl(lsta,:).TypeOrigin,tbl(lsta,:).TypeDest));
                            lsta(lst2)=[];
                        else
                            lst2=find(tbl(lst,:).TypeOrigin==tbl(lst,:).TypeDest);
                            lst(lst2)=[];

                            lst2=find(tbl(lstc(lsta),:).TypeOrigin==tbl(lstc(lsta),:).TypeDest);
                            lsta(lst2)=[];

                        end

                    end
                    data(i).probe.connections=data(i).probe.connections(lst,:);
                    data(i).R=data(i).R(lst);
                    if(~isempty(data(i).ZstdErr))
                        data(i).ZstdErr=data(i).ZstdErr(lsta,:,:);
                    end
                    
                    
                else
                    switch class(data)
                        case {'nirs.core.Data'}
                            tdata = data(i).probe.link.type;

                        case {'nirs.core.ChannelStats','nirs.core.ChannelFStats','nirs.core.ImageStats'}
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

                        case {'nirs.core.ChannelFStats'}
                            data(i).variables = data(i).variables(keep_type,:);
                            data(i).F = data(i).F(keep_type);
                        case {'nirs.core.ImageStats'}
                            data(i).variables = data(i).variables(keep_type,:);
                            data(i).beta = data(i).beta(keep_type);
                            data(i).covb_chol = data(i).covb_chol(keep_type,:);
                            data(i).typeII_StdE = data(i).typeII_StdE(keep_type);
                        case {'nirs.core.sFCStats'}
                            data(i).R = data(i).R(keep_type,keep_type,:);
                            if(~isempty(data(i).ZstdErr))
                                data(i).ZstdErr = data(i).ZstdErr(keep_type,keep_type,:,:);
                            end
                    end
                end
            end
        end
    end
end
