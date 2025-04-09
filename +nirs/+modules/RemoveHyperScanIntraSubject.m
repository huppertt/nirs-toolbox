classdef RemoveHyperScanIntraSubject < nirs.modules.AbstractModule
    %% RemoveHyperScanIntraSubject - Removes all intra subject (within) from hyperscaning
    
    

    methods
        function obj = RemoveHyperScanIntraSubject( prevJob )
            if nargin > 0; obj.prevJob = prevJob; end
            obj.name = 'Remove within subject connections';
        end
        
        function data = runThis( obj, data )
            
            
            for i = 1:numel(data)
                if(isa(data(i),'nirs.core.sFCStats'))
                    if(isa(data(i).probe,'nirs.core.ProbeHyperscan') || ...
                            isa(data(i).probe,'nirs.core.ProbeHyperscan1020'))
                        ll=data(i).probe.link;
                        lst=find(strcmp(ll.SubjectLabelOrigin,ll.SubjectLabelDest));
                        data(i).R(lst)=[];
                        data(i).probe.connections(lst,:)=[];
                        if(~isempty(data(i).ZstdErr))
                            data(i).ZstdErr(lst,:,:)=[];
                        end
                    end
                end
            end
        end

    end
end
