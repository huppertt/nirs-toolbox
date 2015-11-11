classdef MixedEffectsConnectivity < nirs.modules.AbstractModule
    %% MixedEffect - Performs group level mixed effects analysis.
    %
    % Options:
    %     formula     - string specifiying regression formula (see Wilkinson notation)
    %     dummyCoding - dummyCoding format for categorical variables (full, reference, effects)
    %     centerVars  - (true or false) flag for whether or not to center numerical variables
    %
    % Example Formula:
    %     % this will calculate the group average for each condition
    %     j = nirs.modules.MixedEffects();
    %     j.formula = 'beta ~ -1 + group:cond + (1|subject)';
    %     j.dummyCoding = 'full';
    
    properties
        formula = 'F ~ -1 + channel';
        dummyCoding = 'full';
        centerVars = true;
    end
    
    methods
        function obj = MixedEffects( prevJob )
            obj.name = 'Mixed Effects Model for Connectivity';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function G = runThis( obj, S )
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            
            if(isempty(demo))
                for idx=1:length(S)
                    S(idx).demographics('fileIdx')=idx;
                end
                demo = nirs.createDemographicsTable( S );
            end
            
            % center numeric variables
            if obj.centerVars
                n = demo.Properties.VariableNames;
                for i = 1:length(n)
                    if all( isnumeric( demo.(n{i}) ) )
                        demo.(n{i}) = demo.(n{i}) - mean( demo.(n{i}) );
                    end
                end
            end
            
            % Let's do this per channel for now
            n=height(S(1).table);
            
              if(strcmp(S(1).type,'Grangers'))
                  fld='Grangers';
              else
                  fld='Z';
              end
            
            D=zeros(length(S),n);
            for i=1:length(S)
                D(i,:)=real(S(i).(fld)(:));
            end
            D(D==Inf)=1/eps(1);
            D(D==-Inf)=-1/eps(1);
            
            formula=obj.formula;
            formula=['corr ' formula(strfind(formula,'~'):end)];
            nRE=length(strfind(obj.formula,'|'));
            
            for idx=1:n
                corr=D(:,idx);
                vars=[demo table(corr)];
                if(nRE>0)
                    lm = fitlme(vars,formula, 'dummyVarCoding',obj.dummyCoding,...
                        'FitMethod', 'ML');
                else
                    lm = fitlm(vars,formula, 'dummyVarCoding',obj.dummyCoding);    
                end
            end
            
%             
%             vars = table();
%             for i = 1:length(S)
%                 tbl=S(i).table;
%                 tbl=[tbl repmat(demo(i,:),height(tbl),1)];
%                 %tbl=sortrows(tbl,{'TypeOrigin','TypeDest'});
%                 LabelsOrig=strcat(repmat('src',height(tbl),1),num2str(tbl.SourceOrigin),...
%                     repmat('det',height(tbl),1), num2str(tbl.DetectorOrigin),tbl.TypeOrigin);
%                 
%                 LabelsDest=strcat(repmat('src',height(tbl),1),num2str(tbl.SourceDest),...
%                     repmat('det',height(tbl),1), num2str(tbl.DetectorDest),tbl.TypeDest);
%                 channel=strcat(LabelsOrig,repmat('_',height(tbl),1),LabelsDest);
%                 vars=[vars; [table(channel) tbl]];
%                 
%             end
%             
%             nRE=length(strfind(obj.formula,'|'));
%             formula=obj.formula;
%             
%             if(strcmp(S(1).type,'Grangers'))
%                 
%                 formula=['Grangers ' formula(strfind(formula,'~'):end)];
%                 lm = fitlme(vars,formula, 'dummyVarCoding',obj.dummyCoding,...
%                     'FitMethod', 'ML', 'CovariancePattern', repmat({'Diagonal'},nRE,1),...
%                     'Verbose',true);
%             else
%                 
%                 vars.Z(vars.Z==Inf)=1/eps(1);
%                 vars.Z(vars.Z==-Inf)=-1/eps(1);
%                 formula=['Z ' formula(strfind(formula,'~'):end)];
%                 lm = fitlme(vars,formula, 'dummyVarCoding',obj.dummyCoding,...
%                     'FitMethod', 'ML', 'CovariancePattern', repmat({'Diagonal'},nRE,1),...
%                     'Verbose',true);
%             end
%             
            %Now sort back out
            G = nirs.core.ConnectivityStats();
            G.description = 'Group Level Connectivity';
            G.type=S(1).type;
            G.probe=S(1).probe;
            
            %  Labels=strcat(repmat('Labels_',length(Labels),1),Labels);
            if(ismember('conditions',vars.Properties.VariableNames))
                nConds=length(unique(vars.conditions));
            else
                nConds=1;
            end
            if(strcmp(G.type,'Grangers'))
                [n,m]=size(S(1).Grangers);
                % lst=find(ismember(lmG.CoefficientNames,Labels));
                Gr=reshape(lm.Coefficients.Estimate,n,m,nConds);
                dfe1=S(1).dfe1; for idx=2:length(S); dfe1=dfe1+S(idx).dfe1; end;
                dfe2=S(1).dfe2; for idx=2:length(S); dfe2=dfe2+S(idx).dfe2; end;
                G.dfe2=dfe2/length(S);
                G.dfe1=dfe1/length(S);
                G=G.GtoF(Gr);
            else
                [n,m]=size(S(1).Pearsons);
                Z=reshape(lm.Coefficients.Estimate,n,m);
                G.Pearsons=tanh(Z);
                dfe2=S(1).dfe2; for idx=2:length(S); dfe2=dfe2+S(idx).dfe2; end;
                G.dfe2=dfe2/length(S);
                
            end
            
        end
    end
    
end