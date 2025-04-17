classdef ConnectivitySymmetry < nirs.modules.AbstractModule
    %% RemoveHyperScanIntraSubject - Removes all intra subject (within) from hyperscaning
    
    
    properties
        force = false; 
    end
    methods
        function obj = ConnectivitySymmetry( prevJob )
            if nargin > 0; obj.prevJob = prevJob; end
            obj.name = 'Prune the connecitvity data to enforce symmetry';
        end

        function data = runThis( obj, data )


            for i = 1:numel(data)
                if(isa(data(i),'nirs.core.sFCStats'))
                    if(isa(data(i).probe,'nirs.core.ProbeHyperscan'))

                        tbl=data(i).table;
                        
                        % Remove all the redundant inter-subject
                        % connections
                        for id=1:height(tbl);
                            org=[tbl.condition{id}, '_', tbl.TypeOrigin{id}, 'S', num2str(tbl.SourceOrigin(id)),...
                                ':D', num2str(tbl.DetectorOrigin(id))];
                            dest=[tbl.condition{id}, '_' tbl.TypeDest{id}, 'S', num2str(tbl.SourceDest(id)),...
                                ':D', num2str(tbl.DetectorDest(id))];
                            if(strcmp(tbl.SubjectLabelOrigin{id},tbl.SubjectLabelDest{id}))
                                str=sort({org dest});
                                vv{id,1}=[tbl.SubjectLabelOrigin{id} '_' str{1} '_' str{2}];
                            else
                                if(data(i).ishypersymm)
                                    str=sort({tbl.SubjectLabelOrigin{id} tbl.SubjectLabelDest{id}});
                                    vv{id,1}=[str{1} '_' str{2} '_' org '_' dest];
                                else
                                    str={org dest};
                                    vv{id,1}=[str{1} '_' str{2} '_' org '_' dest];
                                end
                            end
                        end;
                        [uvv,idx]=unique(vv);

                        for id=1:length(idx);
                            lst=find(ismember(vv,uvv{id}));
                            maxV(id)=var(atanh(tbl.R(lst)));
                            R(id,1)=tanh(mean(atanh(tbl.R(lst))));
                        end;

                        data(i).probe.connections=data(i).probe.connections(idx,:);
                        data(i).R=data(i).R(idx);


                        if(max(maxV)>eps(1000));
                            if(obj.force)
                                warning('Connectivity matrix is not symetric!');
                                data(i).R=R;
                            else
                                error('Connectivity matrix is not symetric!');
                            end
                        end


                        if(~isempty(data(i).ZstdErr))
                            data(i).ZstdErr=data(i).ZstdErr(idx,:);
                        end
                        
                       


                    else
                        tbl=data(i).table;
                        for id=1:height(tbl);
                            org=[tbl.condition{id}, '_', tbl.TypeOrigin{id}, 'S', num2str(tbl.SourceOrigin(id)),...
                                ':D', num2str(tbl.DetectorOrigin(id))];
                            dest=[tbl.condition{id}, '_' tbl.TypeDest{id}, 'S', num2str(tbl.SourceDest(id)),...
                                ':D', num2str(tbl.DetectorDest(id))];

                          str=sort({org dest});
                            vv{id,1}=[str{1} '_' str{2}];
                        end;

                        [uvv,idx]=unique(vv);

                        for id=1:length(idx);
                            lst=find(ismember(vv,uvv{id}));
                            maxV(id)=var(tbl.R(lst));
                            R(id,1)=mean(tbl.R(lst));
                        end;

                        data(i).probe.connections=data(i).probe.connections(idx,:);
                        data(i).R=data(i).R(idx);


                        if(max(maxV)>eps(1000));
                            if(obj.force)
                                warning('Connectivity matrix is not symetric!');
                                data(i).R=R;
                            else
                                error('Connectivity matrix is not symetric!');
                            end
                        end


                        if(~isempty(data(i).ZstdErr))
                            data(i).ZstdErr=data(i).ZstdErr(idx,:);
                        end
                    end
                end
            end
        end

    end
end
