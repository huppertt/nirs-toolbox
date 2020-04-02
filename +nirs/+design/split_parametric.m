function d = split_parametric(d,formula,centervariables)
% This function splits a parametric stimulus design creating stimEvents
% variables (and placing them into the data) according to the formula. 
% The default formula is '<Cond>*(1+Amp)' which means a stimulus with 
% parametric modulation on amplitude will be split to a constant amplitude 
% regressor and a linear term.
%  Other options include quadratic '<Cond>*(1+Amp+Amp^2)'
%  or time varing '<Cond>*(1+time)' models
% <cond> needs to be the name of the condition (use ? for wildcard)
% e.g. "SOT?" means all stim starting with SOT whereas ? means all stim 


if(nargin<2)
    formula = '?*(1+amp)';
end

if(nargin<3)
    centervariables=true;
end

if(length(d)>1)
    for idx=1:length(d)
        d(idx)=nirs.design.split_parametric(d(idx),formula,centervariables);
    end
    return
end

flag=1*centervariables;
    formula(isspace(formula))=[];

cond = formula(1:min(strfind(formula,'*'))-1);
if(isempty(cond))
    cond = formula;
end

eqn =  formula(min(strfind(formula,'*'))+1:end);

cond(strfind(cond,'<'))=[];
cond(strfind(cond,'>'))=[];

% Start by dealing with the wildcard/multiple condition versions of the
% formula
if(strcmp(cond,'?'))
    conds=d.stimulus.keys;
    for idx=1:length(conds)
        d=nirs.design.split_parametric(d,[conds{idx} '*' eqn],centervariables);
    end
    return;
end

if(~isempty(strfind(cond,'?')))
    conds=d.stimulus.keys;
    cond(strfind(cond,'?'))=[];
    for idx=1:length(conds)
        if(~isempty(strfind(cond,conds{idx})))
            d=nirs.design.split_parametric(d(idx),[conds{idx} '*' eqn],centervariables);
        end
    end
    return;
end

cond(strfind(cond,'<'))=[];
cond(strfind(cond,'>'))=[];



% Its probably just easier to brute force our way through this.  I hate
% coding like this, but ...

stim = d.stimulus(cond);
if(isa(stim,'nirs.design.StimulusEvents'))
    
    stim.name=cond;
    newstim={};
    while(1)
        cnt=0;
        
               
        
        if(isempty(eqn))
            newstim{end+1}=stim;
            newstim{end}.name = stim.name;
            newstim{end}.amp = ones(size(stim.amp));
            break;
        end
        
        if(~isempty(strfind(eqn,'1')))
            newstim{end+1}=stim;
            newstim{end}.name = stim.name;
            newstim{end}.amp = ones(size(stim.amp));
            eqn(strfind(eqn,'1'))=[];
            cnt=cnt+1;
        end
        
        if(~isempty(strfind(lower(eqn),'+amp^')))
            exp = str2num(eqn(strfind(lower(eqn),'amp^')+4));
            newstim{end+1}=stim;
            newstim{end}.name = [stim.name ':amp^' num2str(exp)];
            newstim{end}.amp = (stim.amp.^exp-flag*mean(stim.amp.^exp));
            newstim{end}.amp=newstim{end}.amp./max(abs(newstim{end}.amp));
            eqn(strfind(lower(eqn),'amp^')+[0:4])=[];
            cnt=cnt+1;
        elseif(~isempty(strfind(lower(eqn),'-amp^')))
            exp = str2num(eqn(strfind(lower(eqn),'amp^')+4));
            newstim{end+1}=stim;
            newstim{end}.name = [stim.name ':amp^' num2str(exp)];
            newstim{end}.amp = (stim.amp.^exp-flag*mean(stim.amp.^exp));
            newstim{end}.amp=-newstim{end}.amp./max(abs(newstim{end}.amp));
            eqn(strfind(lower(eqn),'amp^')+[0:4])=[];
            cnt=cnt+1;
        elseif(~isempty(strfind(lower(eqn),'+amp')))
            newstim{end+1}=stim;
            newstim{end}.name = [stim.name ':amp'];
            newstim{end}.amp = (stim.amp-flag*mean(stim.amp));
            newstim{end}.amp=newstim{end}.amp./max(abs(newstim{end}.amp));
            eqn(strfind(lower(eqn),'amp')+[0:2])=[];
            cnt=cnt+1;
                    elseif(~isempty(strfind(lower(eqn),'-amp')))
            newstim{end+1}=stim;
            newstim{end}.name = [stim.name ':amp'];
            newstim{end}.amp = (stim.amp-flag*mean(stim.amp));
            newstim{end}.amp=-newstim{end}.amp./max(abs(newstim{end}.amp));
            eqn(strfind(lower(eqn),'amp')+[0:2])=[];
            cnt=cnt+1;
        end
        
        if(~isempty(strfind(lower(eqn),'+dur^')))
            exp = str2num(eqn(strfind(lower(eqn),'dur^')+4));
            newstim{end+1}=stim;
            newstim{end}.name = [stim.name ':dur^' num2str(exp)];
            newstim{end}.amp = ones(size(stim.amp)).*(stim.dur.^exp-flag*mean(stim.dur.^exp));
            newstim{end}.amp=newstim{end}.amp./max(abs(newstim{end}.amp));
            eqn(strfind(lower(eqn),'dur')+[0:4])=[];
            cnt=cnt+1;
        elseif(~isempty(strfind(lower(eqn),'-dur^')))
            exp = str2num(eqn(strfind(lower(eqn),'dur^')+4));
            newstim{end+1}=stim;
            newstim{end}.name = [stim.name ':dur^' num2str(exp)];
            newstim{end}.amp = ones(size(stim.amp)).*(stim.dur.^exp-flag*mean(stim.dur.^exp));
            newstim{end}.amp=-newstim{end}.amp./max(abs(newstim{end}.amp));
            eqn(strfind(lower(eqn),'dur')+[0:4])=[];
            cnt=cnt+1;
        elseif(~isempty(strfind(lower(eqn),'+dur')))
            newstim{end+1}=stim;
            newstim{end}.name = [stim.name ':dur'];
            newstim{end}.amp = ones(size(stim.amp)).*(stim.dur-flag*mean(stim.dur));
            newstim{end}.amp=newstim{end}.amp./max(abs(newstim{end}.amp));
            eqn(strfind(lower(eqn),'dur')+[0:2])=[];
            cnt=cnt+1;
             elseif(~isempty(strfind(lower(eqn),'-dur')))
            newstim{end+1}=stim;
            newstim{end}.name = [stim.name ':dur'];
            newstim{end}.amp = ones(size(stim.amp)).*(stim.dur-flag*mean(stim.dur));
            newstim{end}.amp=-newstim{end}.amp./max(abs(newstim{end}.amp));
            eqn(strfind(lower(eqn),'dur')+[0:2])=[];
            cnt=cnt+1;
        end
        
        
        if(~isempty(strfind(lower(eqn),'+time^')))
            exp = str2num(eqn(strfind(lower(eqn),'time^')+5));
            newstim{end+1}=stim;
            newstim{end}.name = [stim.name ':time^' num2str(exp)];
            newstim{end}.amp = ones(size(stim.amp)).*(stim.onset.^exp-flag*mean(stim.onset.^exp));
            newstim{end}.amp=sign*newstim{end}.amp./max(abs(newstim{end}.amp));
            eqn(strfind(lower(eqn),'time^')+[0:5])=[];
            cnt=cnt+1;
        elseif(~isempty(strfind(lower(eqn),'-time^')))
            exp = str2num(eqn(strfind(lower(eqn),'time^')+5));
            newstim{end+1}=stim;
            newstim{end}.name = [stim.name ':time^' num2str(exp)];
            newstim{end}.amp = ones(size(stim.amp)).*(stim.onset.^exp-flag*mean(stim.onset.^exp));
            newstim{end}.amp=-newstim{end}.amp./max(abs(newstim{end}.amp));
            eqn(strfind(lower(eqn),'time^')+[0:5])=[];
            cnt=cnt+1;
        elseif(~isempty(strfind(lower(eqn),'+time')))
            newstim{end+1}=stim;
            newstim{end}.name = [stim.name ':time'];
            newstim{end}.amp = ones(size(stim.amp)).*(stim.onset-flag*mean(stim.onset));
            newstim{end}.amp=newstim{end}.amp./max(abs(newstim{end}.amp));
            eqn(strfind(lower(eqn),'time')+[0:3])=[];
            cnt=cnt+1;
            elseif(~isempty(strfind(lower(eqn),'-time')))
            newstim{end+1}=stim;
            newstim{end}.name = [stim.name ':time'];
            newstim{end}.amp = ones(size(stim.amp)).*(stim.onset-flag*mean(stim.onset));
            newstim{end}.amp=-newstim{end}.amp./max(abs(newstim{end}.amp));
            eqn(strfind(lower(eqn),'time')+[0:3])=[];
            cnt=cnt+1;
        end
        
        if(cnt==0)
            break;
        end
        if(isempty(eqn))
            break
        end
        
    end
    
    for idx=1:length(newstim)
        d.stimulus(newstim{idx}.name)=newstim{idx};
    end
end

return;





