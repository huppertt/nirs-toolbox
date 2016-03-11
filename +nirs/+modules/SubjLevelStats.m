classdef SubjLevelStats < nirs.modules.AbstractModule
    %% SubjLevelStats - Performs subject level analysis to combine all files for each subject into one stats object.
    %
    % Options:
    %
    properties
        sortfield='subject';
        
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
           
            if(~ismember(demo.Properties.VariableNames,obj.sortfield))
                warning('Sort field not contained in demograophics');
                G=S;
                return
            end
            
            subjects = unique(demo.(obj.sortfield));
            
            if(~iscell(subjects)); for i=1:length(subjects); ss{i}=subjects(i); end; subjects=ss; end;
            
            if(length(subjects)>1)
                G = nirs.core.ChannelStats();
                for i=1:length(subjects)
                    lst=find(ismember(demo.(obj.sortfield),subjects{i}));
                    G(i)=obj.runThis(S(lst));
                end
                return
            end
                       
            % preallocate group stats
            G = nirs.core.ChannelStats();
            
            %% loop through files
            W = sparse([]);
            iW = sparse([]);
            
            b = [];
            vars = table();
            C=[];
            for i = 1:length(S)
                % coefs
                b = [b; S(i).beta];
                
                C=blkdiag(C,S(i).covb);
                
                % table of variables
                file_idx = repmat(i, [size(S(i).beta,1) 1]);
                
                if(~isempty(demo))
                    vars = [vars;
                        [table(file_idx) S(i).variables repmat(demo(i,:), [size(S(i).beta,1) 1])]
                        ];
                else
                    vars = [vars; ...
                        [table(file_idx) S(i).variables]];
                end
            end
            
            % sort
            [vars, idx] = sortrows(vars, {'source', 'detector', 'type'});
            
            % list for first source
            [sd, ~,lst] = unique(table(vars.source, vars.detector, vars.type), 'rows', 'stable');
            sd.Properties.VariableNames = {'source', 'detector', 'type'};
            
            %% design mats
            tmp = vars(lst == 1, :);
            
            beta = randn(size(tmp,1), 1);
            
            warning('off','stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_PerfectFit');
            lm1=fitlme([table(beta) tmp], 'beta ~ -1 + cond','dummyVarCoding','full');
                        
            X = lm1.designMatrix('Fixed');
            nchan = max(lst);
            X = kron(speye(nchan), X);
            
            %% put them back in the original order
            vars(idx,:) = vars;
            X(idx, :)   = X;
            beta        = b; % already in correct order
            
            X=X./(ones(size(X,1),1)*sum(X,1));
            
            cnames = lm1.CoefficientNames(:);
            for idx=1:length(cnames);
                cnames{idx}=cnames{idx}(max([0 min(strfind(cnames{idx},'_'))])+1:end);
                %if(cnames{idx}(1)=='_'); cnames{idx}(1)=[]; end;
            end;
            cnames = repmat(cnames, [nchan 1]);
            
            G.beta = X' * beta;
            G.covb = X'*C*X;
            G.dfe = sum(vertcat(S.dfe));
            G.probe=S(1).probe;
            
           
            
            sd = repmat(sd, [length(unique(cnames)) 1]);
            sd = sortrows(sd, {'source', 'detector', 'type'});
            
            G.variables = [sd table(cnames)];
            G.variables.Properties.VariableNames{4} = 'cond';
            G.description = 'Subject level Model';
            
            G.demographics=S(1).demographics;
            
            
            
                        
        end
    end
    
end
