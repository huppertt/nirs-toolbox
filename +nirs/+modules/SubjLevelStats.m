classdef SubjLevelStats < nirs.modules.AbstractModule
    %% SubjLevelStats - Performs subject level analysis to combine all files for each subject into one stats object.
    %
    % Options:
    %
    properties
        sortfield='subject';
        formula='beta ~ -1 + cond';
        robust=false;
        verbose=false;
        include_diagnostics=false;
    end
    
    
    methods
        function obj = SubjLevelStats( prevJob )
            obj.name = 'Subject level Effects Model';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function G = runThis( obj, S )
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
           
            if(~isempty(obj.sortfield))
                if(~ismember(demo.Properties.VariableNames,obj.sortfield))
                    warning('Sort field not contained in demograophics');
                    G=S;
                    return
                end
                
                subjects = demo.(obj.sortfield);
                if(ischar(subjects)); subjects=cellstr(subjects); end;
                subjects = unique(subjects);
                
                if(length(subjects)>1)
                    G = nirs.core.ChannelStats();
                    for i=1:length(subjects)
                        lst=find(ismember(demo.(obj.sortfield),subjects{i}));
                        
                        cond=unique(nirs.getStimNames(S(lst)));
                        hascond = false(length(lst),length(cond));
                        for j=1:length(lst)
                            hascond(j,:) = ismember(cond,nirs.getStimNames(S(lst(j))));
                        end
                        lstSingleCond = find(sum(hascond*1,1)==1);
                        j=nirs.modules.KeepStims;
                        j.listOfStims={cond{lstSingleCond}};
                        SS = j.run(S(lst));
                        
                        j=nirs.modules.DiscardStims;
                        j.listOfStims={cond{lstSingleCond}};
                        j=nirs.modules.RemoveStimless(j);
                        Scommon = j.run(S(lst));
                        
                        if(~isempty( Scommon))
                            disp([ 'running subject ' subjects{i}]);
                                       G1=obj.runThis( Scommon);
                            if(~isempty(SS))
                                SS(end+1)=G1;
                            else
                                SS=G1;
                            end
                        end
                        G(i,1)=nirs.math.combineStats(SS);
                        
                    end
                    disp(' ');
                    return
                end
            end
                    
            % Return if there is nothing to average
            if(height(demo)==1)
                G=S;
                return;
            end
            
       
            
            j=nirs.modules.MixedEffects;
            j.formula=obj.formula;
            j.verbose=obj.verbose;
            j.include_diagnostics=obj.include_diagnostics;
            j.robust=obj.robust;
            
            G = j.run(S);
            
%             
%             % preallocate group stats
%             G = nirs.core.ChannelStats();
%             
%             %% loop through files
%             W = sparse([]);
%             iW = sparse([]);
%             
%             b = [];
%             vars = table();
%             C=[];
%             for i = 1:length(S)
%                 % coefs
%                 b = [b; S(i).beta];
%                 
%                 C=blkdiag(C,S(i).covb);
%                 
%                 % table of variables
%                 file_idx = repmat(i, [size(S(i).beta,1) 1]);
%                 
%                 if(~isempty(demo))
%                     vars = [vars;
%                         [table(file_idx) S(i).variables repmat(demo(i,:), [size(S(i).beta,1) 1])]
%                         ];
%                 else
%                     vars = [vars; ...
%                         [table(file_idx) S(i).variables]];
%                 end
%             end
%             
%             % sort
%             [vars, idx] = nirs.util.sortrows(vars, {'source', 'detector', 'type'});
%             
%             % list for first source
%             [sd, ~,lst] = nirs.util.uniquerows(table(vars.source, vars.detector, vars.type));
%             sd.Properties.VariableNames = {'source', 'detector', 'type'};
%             
%             %% design mats
%             tmp = vars(lst == 1, :);
%             
%             beta = randn(size(tmp,1), 1);
%             
%             warning('off','stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_PerfectFit');
%             lm1=fitlme([table(beta) tmp], obj.formula,'dummyVarCoding','full');
%                         
%             X = lm1.designMatrix('Fixed');
%             nchan = max(lst);
%             X = kron(speye(nchan), X);
%             
%             %% put them back in the original order
%             vars(idx,:) = vars;
%             X(idx, :)   = X;
%             beta        = b; % already in correct order
%             
%             X=X./(ones(size(X,1),1)*sum(X,1));
%             
%             cnames = lm1.CoefficientNames(:);
%             for idx=1:length(cnames);
%                 cnames{idx}=cnames{idx}(max([0 min(strfind(cnames{idx},'_'))])+1:end);
%                 %if(cnames{idx}(1)=='_'); cnames{idx}(1)=[]; end;
%             end;
%             cnames = repmat(cnames, [nchan 1]);
%             
%             G.beta = X' * beta;
%             G.covb = X'*C*X;
%             G.dfe = sum(vertcat(S.dfe));
%             G.probe=S(1).probe;
%             
%            
%             
%             sd = repmat(sd, [length(unique(cnames)) 1]);
%             sd = nirs.util.sortrows(sd, {'source', 'detector', 'type'});
%             
%             G.variables = [sd table(cnames)];
%             G.variables.Properties.VariableNames{4} = 'cond';
%             G.description = 'Subject level Model';
%             
%             G.demographics=S(1).demographics;
%             
%             n={}; b={}; cnt=1;
%             for i=1:length(S)
%                 for j=1:S(i).basis.stim.count;
%                     n{cnt}=S(i).basis.stim.values{j}.name;
%                     b{cnt}=S(i).basis.stim.values{j};
%                     cnt=cnt+1;
%                 end
%             end
%             [~,j]=unique(n);
%             G.basis=S(1).basis;
%             G.basis.stim=Dictionary;
%             for i=1:length(j)
%                 G.basis.stim(n{j(i)})=b{j(i)};
%             end
            
            
                        
        end
    end
    
end
