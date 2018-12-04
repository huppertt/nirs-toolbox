classdef AnovaN < nirs.modules.AbstractModule
    %N-way ANOVA module based on the Matlab anovan function.
  
    properties
        depvar ='beta';
        variables={'cond','subject','group'};
        model = 'linear';
        sstype=3;
    end

    methods
        function obj = AnovaN( prevJob )
            obj.name = 'Anova Model';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function G = runThis( obj, S )
            
              
            % demographics info
            demo = nirs.createDemographicsTable( S );
            vars = table(); b=[];
            for i = 1:length(S)
                if(strcmp(obj.depvar,'tstat'))
                    b = [b; S(i).tstat];
                else
                    b = [b; S(i).beta];
                end
                file = repmat(cellstr(num2str(i)), [size(S(i).beta,1) 1]);
                 if(~isempty(demo))
                    vars = [vars;
                        [table(file) S(i).variables repmat(demo(i,:), [size(S(i).beta,1) 1])]
                        ];
                else
                    vars = [vars; ...
                        [table(file) S(i).variables]];
                 end
            end
           % list for first source
           if(~ismember('source',vars.Properties.VariableNames) & ...
                    ismember('ROI',vars.Properties.VariableNames))
                [sd, ~,lst] = nirs.util.uniquerows(table(vars.ROI, vars.type));
            sd.Properties.VariableNames = {'ROI', 'type'};   
           else
            [sd, ~,lst] = nirs.util.uniquerows(table(vars.source, vars.detector, vars.type));
            sd.Properties.VariableNames = {'source', 'detector', 'type'};    
           end
           
            vars=vars(:,find(ismember(vars.Properties.VariableNames,obj.variables)));
            
             % preallocate group stats
            G = nirs.core.ChannelFStats();
            
            for i=1:height(sd)
                ll=find(lst==i);
                
                c=false(size(vars,2),1);  %catagorical or not
                for j=1:size(vars,2); 
                    L{j}=table2cell(vars(ll,j)); 
                    try;
                        if(isnumeric(vertcat(L{j}{:})))
                            L{j}=vertcat(L{j}{:});
                            c(j)=true;
                        end
                    end
                    
                end;
               c=find(c);
                
               [p,atab,stats,terms]=anovan(b(ll),L,'sstype',obj.sstype,...
                   'model',obj.model,'display','off','varnames',vars.Properties.VariableNames,'continuous',c);
            
                for j=1:size(atab,1)-3; 
                    F(i,j)=atab{1+j,6};
                    df1(i,j)=atab{1+j,3};
                    df2(i,j)=length(ll)-df1(i,j)-1;
                    name{j}=atab{1+j,1};
                end
            end
            
            G.variables=table;
            for j=1:length(name)
                if(all(F(:,j)==0))
                    continue;
                end
                 G.variables=[ G.variables; sd table(repmat({name{j}},height(sd),1),'VariableName',{'cond'})];
                 G.F=[G.F; F(:,j)];
                 G.df1=[G.df1; df1(:,j)];
                 G.df2=[G.df2; df2(:,j)];
                 
            end
            G.probe=S(1).probe;
            G.probe.link=sd;
            G = G.sorted();
            G.description = 'nway-ANOVA Model';
            
            G.demographics = nirs.util.combine_demographics(...
                nirs.createDemographicsTable(S));
            
        end
            
    end
    
end

