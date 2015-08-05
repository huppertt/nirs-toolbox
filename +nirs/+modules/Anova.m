classdef Anova < nirs.modules.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        formula = 'beta ~ 1 + group*cond + (1|subject)';
        centerVars = false;
    end

    methods
        function obj = Anova( prevJob )
            obj.name = 'Anova Model';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function G = runThis( obj, S )
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            
            if obj.centerVars
                % center numeric variables
                n = demo.Properties.VariableNames;
                for i = 1:length(n)
                    if all( isnumeric( demo.(n{i}) ) )
                        demo.(n{i}) = demo.(n{i}) - mean( demo.(n{i}) );
                    end
                end
            end
            
            % preallocate group stats
            G = nirs.core.ChannelFStats();
            
            %% loop through files
            w = [];
            b = [];
            vars = table();
            for i = 1:length(S)
                % coefs
                b = [b; S(i).beta];
                
                % weights
                w = [w; 1./sqrt(diag(S(i).covb))];
                
                % table of variables
                file_idx = repmat(i, [size(S(i).beta,1) 1]);
                
                vars = [vars;
                    [table(file_idx) S(i).variables repmat(demo(i,:), [size(S(i).beta,1) 1])]
                    ];
            end
            
            % sort
            [vars, idx] = sortrows(vars, {'source', 'detector', 'type', 'cond'});
            
            % list for first source
            [sd, ~,lst] = unique(table(vars.source, vars.detector, vars.type), 'rows', 'stable');
            sd.Properties.VariableNames = {'source', 'detector', 'type'};
            
            %% loop through
            variables = table([],[],[],[], 'VariableNames', {'source', 'detector', 'type', 'cond'});
            F = []; df1 = []; df2 = [];
            for iChan = 1:max(lst)
                
                tmp = vars(lst == iChan, :);

                beta = b(lst == iChan);
                
                lm = fitlme([table(beta) tmp], obj.formula, 'dummyVarCoding',...
                        'reference', 'FitMethod', 'ML', 'CovariancePattern', 'Isotropic', ...
                        'Weights', w(lst==iChan));
                    
                a = lm.anova();
                
                F = [F; a.FStat];
                df1 = [df1; a.DF1];
                df2 = [df2; a.DF2];
                
                cond = a.Term;
                variables = [variables; 
                    [repmat(sd(iChan,:), [length(a.FStat) 1]) table(cond)]];
            end
            
            G.F = F;
            G.df1 = df1;
            G.df2 = df2;
            G.variables = variables;
            G.probe = S(1).probe;
            
            G = G.sorted();
        end
    end
    
end

