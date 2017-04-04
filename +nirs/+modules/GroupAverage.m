classdef GroupAverage < nirs.modules.AbstractModule
    %% Group Average - Performs group level mixed effects analysis.
    % This function is essentially the same as the Mixed Effects model but
    % does not support the random  effect terms.  
    % This uses the robust fit code to do the average.  If you don't care 
    % about the random effects terms, then this is several
    % fold faster then the Mixed Effects model.   
    % 
    % Options:
    %     formula     - string specifiying regression formula (see Wilkinson notation)
    %     dummyCoding - dummyCoding format for categorical variables (full, reference, effects)
    %     centerVars  - (true or false) flag for whether or not to center numerical variables
    %     robustfit  - (true or false) flag to use the robust fit code  
    %
    % Example Formula:
    %     % this will calculate the group average for each condition
    %     j = nirs.modules.MixedEffects();
    %     j.formula = 'beta ~ -1 + group:cond';
    %     j.dummyCoding = 'full';
    
    properties
        formula = 'beta ~ -1 + group:cond';
        dummyCoding = 'full';
        centerVars = true;
        include_diagnostics=false;
        robustfit=false;
    end
    
    methods
        function obj = GroupAverage( prevJob )
            obj.name = 'Group Average Model';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function G = runThis( obj, S )
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            
            % center numeric variables
            if obj.centerVars
                n = demo.Properties.VariableNames;
                for i = 1:length(n)
                    if all( isnumeric( demo.(n{i}) ) )
                        demo.(n{i}) = demo.(n{i}) - mean( demo.(n{i}) );
                    end
                end
            end
            
            if(isa(S,'nirs.core.ChannelStats'))
            % preallocate group stats
            G = nirs.core.ChannelStats();
            
            %% loop through files
            W = sparse([]);
           % iW = sparse([]);
            dWTW = [];
            
            b = [];
            vars = table();
            for i = 1:length(S)
                % coefs
                disp([num2str(i) ' of ' num2str(length(S))])
                b = [b; S(i).beta];
                
                % whitening transform
                

                    L=chol(S(1).covb,'lower');
                    L=sparse(L);
                    w=inv(L);
%                     [u, s, ~] = svd(S(i).covb, 'econ');
%                     %W = blkdiag(W, diag(1./diag(sqrt(s))) * u');
%                     %W = blkdiag(W, pinv(s).^.5 * u');
%                     w=sparse(diag(sqrt(1./(diag(s)+eps(1))))) * u';
                    
                 W = sparse(blkdiag(W, w));
                dWTW = [dWTW; sqrt(diag(w'*w))];
                %iW = blkdiag(iW, u*sqrt(s) );
                
                
                %                L = chol(S(i).covb,'upper');
                %                W = blkdiag(W,pinv(L));
                
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
            
            probe=S(1).probe;
            %clear S;
            
            % sort
            [vars, idx] = sortrows(vars, {'source', 'detector', 'type'});
            
            % list for first source
            [sd, ~,lst] = unique(table(vars.source, vars.detector, vars.type), 'rows', 'stable');
            sd.Properties.VariableNames = {'source', 'detector', 'type'};
            
            %% design mats
            tmp = vars(lst == 1, :);
            
            beta = randn(size(tmp,1), 1);
            
            nRE=max(1,length(strfind(obj.formula,'|')));
            lm1 = fitlme([table(beta) tmp], obj.formula,'dummyVarCoding',...
                obj.dummyCoding);
            
            X = lm1.designMatrix('Fixed');
            nchan = max(lst);
            
            X = kron(speye(nchan), X);
            
            %% put them back in the original order
            vars(idx,:) = vars;
            X(idx, :)   = X;
            beta        = b; % already in correct order
            
            %% check weights
           % dWTW = sqrt(diag(W'*W));
            m = median(dWTW);
            
            %W(dWTW > 100*m,:) = 0;
            lstBad=find(dWTW > 100*m);
            
            W(lstBad,:)=[];
            W(:,lstBad)=[];
            X(lstBad,:)=[];
            beta(lstBad,:)=[];
            %% Weight the model
                        
            X    = W*X;
            beta = W*beta;
             
            
            
            %% fit the model
            if(obj.robustfit)
                [G.beta,stats]=nirs.math.robustfit(full(X),beta,[],[],'off');
                G.covb=stats.covb;
                G.dfe= stats.dfe;
            else
                nobs = length(beta);
                [beta,X]=qr(X,beta,0);
                xtx=X'*X;
                xtxi=inv(xtx);
                xbeta=X'*beta;
                
                G.beta = xtxi*xbeta;
             
                p = length(G.beta);
                G.dfe = nobs-p;
                
                sse = norm(beta-X*G.beta)^2;
                mse = sse./G.dfe;
               
                G.covb=xtxi*mse;
                
            end
         
            cnames = lm1.CoefficientNames(:);
            for idx=1:length(cnames);
                cnames{idx}=cnames{idx}(max([0 min(strfind(cnames{idx},'_'))])+1:end);
                %if(cnames{idx}(1)=='_'); cnames{idx}(1)=[]; end;
            end;
            cnames = repmat(cnames, [nchan 1]);
            
            %% output
            %G.beta       = lm2.Coefficients.Estimate;
            %G.covb       = lm2.CoefficientCovariance;
           % G.dfe        = lm2.DFE;
            G.probe      = probe;
            
            sd = repmat(sd, [length(unique(cnames)) 1]);
            sd = sortrows(sd, {'source', 'detector', 'type'});
            
            G.variables = [sd table(cnames)];
            G.variables.Properties.VariableNames{4} = 'cond';
            G.description = ['Mixed Effects Model: ' obj.formula];
           
            G.basis=S(1).basis;
            
            
            
            if(obj.include_diagnostics)
                %Create a diagnotistcs table of the adjusted data
                yproj = beta;
                yproj=inv(W)*yproj;
                
                [sd, ~,lst] = unique(table(vars.source, vars.detector, vars.type), 'rows', 'stable');
                
                btest=[];
                for idx=1:max(lst)                   
                        ll=find(lst == idx);
                        tmp = vars(ll, :);
                        beta = yproj(ll);
                        w=full(dWTW(ll));
                        
                        mdl{idx} = fitlm([table(beta) tmp], lm1.Formula.FELinearFormula,'weights',w.^2,...
                            'dummyVarCoding','full');
                       
                        btest=[btest; mdl{idx}.Coefficients.Estimate];
                                
                end
                mdl=reshape(repmat(mdl,[length(unique(cnames)),1]),[],1);
                G.variables.model=mdl;
            
            end
            else
                %Deal with averaging time courses together
                
                % First, let's just create a giant table to do this
                tbl=table;
                demo=nirs.createDemographicsTable(S);
                for i=1:length(S)
                    data={};
                    for j=1:size(S(i).data,2)
                        data{j,1}=S(i).data(:,j);
                    end
                    
                    tbl=[tbl; [repmat(demo(i,:),size(S(i).data,2),1) S(i).probe.link table(data)]];
                end
                for i=1:height(tbl)
                    ntps(i)=length(tbl.data{i});
                end
                ntps=max(ntps);
                for i=1:height(tbl)
                    tbl.data{i}(end+1:ntps)=NaN;
                end
                [link,~,idx]=unique(table(tbl.source,tbl.detector,tbl.type,'VariableNames',{'source','detector','type'}));
                data=zeros(ntps,height(link));
                for i=1:height(link)
                    lst=find(idx==i);
                    for j=1:length(lst)
                        d(:,j)=tbl.data{lst(j)};
                    end
                    e=std(d,[],2)/sqrt(size(d,2));
                    d=mean(d,2);
                    data(:,i)=d+sqrt(-1)*e;
                end
                G=nirs.core.Data;
                G.description='Group Average model from core.Data';
                G.probe=S(1).probe;
                G.probe.link=link;
                G.data=data;
                G.time=[0:ntps-1]/S(1).Fs;
                G.Fm=S(1).Fm;
            end
            
                        
        end
    end
    
end
