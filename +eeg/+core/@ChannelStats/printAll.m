function printAll( obj, vtype, vrange, thresh, folder, ext )
    %% PRINTALL - draws and saves figures to a folder
    % 
    % Args:
    %     vtype, vrange, thresh  -  same as draw function
    %     folder  -  directory to save to
    %     ext     -  file extension of images (eps, tif, or jpg)

    n = length(obj.conditions);
    
    I = eye(n);
    
    if ~exist(folder, 'dir')
       mkdir(folder); 
    end
    
    utypes = unique(obj.variables.type, 'stable');
    if isnumeric(utypes)
        utypes = arrayfun(@(x) {num2str(x)}, utypes);
    end
    
    for i = 1:length(obj.conditions)
        obj.ttest(I(i,:)).draw(vtype, vrange, thresh);
        for j = length(utypes):-1:1
            f=findobj('name',[utypes{j} ' : ' obj.conditions{i}]);
            figure(f(1));
            set(f, 'PaperPositionMode', 'auto')
            
            if strcmp(ext, 'eps')
                ptype = '-depsc';
            elseif strcmp(ext, 'tif') || strcmp(ext, 'tiff')
                ptype = '-dtiff';
            elseif strcmp(ext, 'jpg') || strcmp(ext, 'jpeg')
                ptype = '-djpeg';
            else
                error('File extension not recognized.')
            end
            
            fname = [obj.conditions{i} '_' utypes{j} '.' ext];
            fname = [folder filesep strjoin(strsplit(fname, ':'), '__')];
            print(f,ptype, fname)
            
            close(f);
        end
        
    end
  
end

