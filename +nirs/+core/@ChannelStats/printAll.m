function printAll( obj, vtype, vrange, thresh, folder, ext, conditionNames )
    %% PRINTALL - draws and saves figures to a folder
    % 
    % Args:
    %     vtype, vrange, thresh  -  same as draw function
    %     folder  -  directory to save to
    %     ext     -  file extension of images (eps, tif, or jpg)
    %     conditionNames - cell array of names for each condition to
    %               override automatic name generation based on conditions with
    %                (should match 1xn conditions)

    n = length(obj.conditions);

    
    
    I = eye(n);
    
    if ~exist(folder, 'dir')
       mkdir(folder); 
    end
    
    utypes = unique(obj.variables.type, 'stable');
    if isnumeric(utypes)
        utypes = arrayfun(@(x) {num2str(x)}, utypes);
    end

    if(nargin<7)
        for i = 1:n
            conditionNames{i}=obj.conditions{i};
        end
    elseif(size(conditionNames,1)~=n)
        error('ConditionNames must match the number of conditions');
    end
    
    for i = 1:n
        %conddata = obj.ttest(I(i,:));
        h = obj.draw(vtype, vrange, thresh,[],{obj.conditions{i}});
        for j = length(utypes):-1:1
            f=intersect(h,findobj('name',[utypes{j} ' : ' obj.conditions{i}]));
            set(f, 'PaperPositionMode', 'auto')
            
            if strcmp(ext, 'eps')
                ptype = '-depsc';
            elseif strcmp(ext, 'tif') || strcmp(ext, 'tiff')
                ptype = '-dtiff';
            elseif strcmp(ext, 'jpg') || strcmp(ext, 'jpeg')
                ptype = '-djpeg';
            elseif strcmp(ext, 'png')
                ptype = '-dpng';
            elseif strcmp(ext, 'fig')
                fname = [obj.conditions{i} '_' utypes{j} '.' ext];
                fname = [folder filesep strjoin(strsplit(fname, ':'), '__')];
                savefig(f, fname)
                
                close(f);
                continue
            else
                error('File extension not recognized.')
            end
            
            fname = [conditionNames{i} '_' utypes{j} '.' ext];
            fname = [folder filesep strjoin(strsplit(fname, ':'), '__')];
            print(f,ptype, fname)
            
            close(f);
        end
        
    end
  
end

