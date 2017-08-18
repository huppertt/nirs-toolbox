function formula = verify_formula(ds,formula,allowRE)
% This function verifies a formula for the fitLME or fitLM models

if(nargin<3)
    allowRE=true;
end
formula=strtrim(formula);
formula(strfind(formula,' '))=[];

treestruct = classreg.regr.LinearMixedFormula(formula);

if(~any(strcmp(treestruct.ResponseName,{'beta','tstat'})) & ...
        ~ismember(treestruct.ResponseName,ds.Properties.VariableNames))
    error('formula should use ''beta'' or ''tstat'' as the response');
end

lst=find(~ismember(treestruct.PredictorNames,ds.Properties.VariableNames));
if(~isempty(lst))
    str='Unknown variables found:';
    for i=1:length(lst)
        str=strcat(str,' | ',char(treestruct.PredictorNames(lst(i))));
    end
    warning(str);
    
    for i=1:length(lst)
        formula(strfind(formula,char(treestruct.PredictorNames(lst(i))))+...
            [0:length(char(treestruct.PredictorNames(lst(i))))-1])=' ';
    end
    formula(strfind(formula,' '))=[];
    
    formula(strfind(formula,'+:')+1)=' ';
    formula(strfind(formula,'-:')+1)=' ';
    if(allowRE)
        formula(strfind(formula,'|)')+[-2:1])=' ';
    end
    formula(strfind(formula,' '))=[];
    if(strcmp(formula(end),'+') | strcmp(formula(end),'-'))
        formula(end)=[];
    end
%     treestruct = classreg.regr.LinearMixedFormula(formula);
%     if(~allowRE)
%         formula=char(treestruct.FELinearFormula);
%     else
%         formula=char(treestruct.char);
%     end
    disp(['Changed formula to : ' formula]);
    
end

