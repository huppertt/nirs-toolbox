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
        formula = 'F ~ -1 + Labels';
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
            
            % center numeric variables
            if obj.centerVars
                n = demo.Properties.VariableNames;
                for i = 1:length(n)
                   if all( isnumeric( demo.(n{i}) ) )
                       demo.(n{i}) = demo.(n{i}) - mean( demo.(n{i}) );
                   end
                end
            end
            
            vars = table();
            for i = 1:length(S)
                tbl=S(i).table;
                %tbl=sortrows(tbl,{'TypeOrigin','TypeDest'});
                LabelsOrig=strcat(repmat('src',height(tbl),1),num2str(tbl.SourceOrigin),...
                    repmat('det',height(tbl),1), num2str(tbl.DetectorOrigin),tbl.TypeOrigin);
                
                LabelsDest=strcat(repmat('src',height(tbl),1),num2str(tbl.SourceDest),...
                    repmat('det',height(tbl),1), num2str(tbl.DetectorDest),tbl.TypeDest);
                Labels=strcat(LabelsOrig,repmat('_',height(tbl),1),LabelsDest);
                vars=[vars; [table(Labels) tbl]];
                
            end
            
            nRE=length(strfind(obj.formula,'|'));
            formula=obj.formula;
            formula=['Grangers ' formula(strfind(formula,'~'):end)];
            lmG = fitlme(vars,formula, 'dummyVarCoding',obj.dummyCoding,...
                'FitMethod', 'ML', 'CovariancePattern', repmat({'diagional'},nRE,1),...
                'Verbose',true);
            
            formula=obj.formula;
            vars.FisherZ(vars.FisherZ==Inf)=1/eps(1);
            vars.FisherZ(vars.FisherZ==-Inf)=-1/eps(1);
            formula=['FisherZ ' formula(strfind(formula,'~'):end)];
            lmZ = fitlme(vars,formula, 'dummyVarCoding',obj.dummyCoding,...
                'FitMethod', 'ML', 'CovariancePattern', repmat({'diagional'},nRE,1),...
                'Verbose',true);
            
            %Now sort back out 
            G = nirs.core.ConnectivityStats();
            G.description = 'Group Level Connectivity';
            G.probe=S(1).probe;
            
            Labels=strcat(repmat('Labels_',length(Labels),1),Labels);
            
            [n,m]=size(S(1).Grangers);
            
            lst=find(ismember(lmG.CoefficientNames,Labels));
            Gr=reshape(lmG.Coefficients.Estimate(lst),n,m);
    
            Z=reshape(lmZ.Coefficients.Estimate,n,m); 
            G.Pearsons=tanh(Z);
            
            dfe1=S(1).dfe1; for idx=2:length(S); dfe1=dfe1+S(idx).dfe1; end;
            dfe2=S(1).dfe2; for idx=2:length(S); dfe2=dfe2+S(idx).dfe2; end;
            G.dfe2=dfe2/length(S);
            G.dfe1=dfe1/length(S);
            G=G.GtoF(Gr);
                
        end
    end
    
end