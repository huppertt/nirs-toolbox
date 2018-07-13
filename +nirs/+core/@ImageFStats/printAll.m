function printAll( obj, fmax, thresh, folder, ext )
    %% PRINTALL - draws and saves figures to a folder
    % 
    % Args:
    %     fmax, thresh  -  same as draw function
    %     folder  -  directory to save to
    %     ext     -  file extension of images (eps, tif, or jpg)

    n = length(obj.conditions);
        
    if ~exist(folder, 'dir')
       mkdir(folder); 
    end
    
    utypes = unique(obj.variables.type, 'stable');
    if isnumeric(utypes)
        utypes = arrayfun(@(x) {num2str(x)}, utypes);
    end
    
    obj.draw(fmax, thresh);
    
    for i = fliplr(1:length(obj.conditions))
        for j = fliplr(1:length(utypes))
            
            set(gcf, 'PaperPositionMode', 'auto')
            
            if strcmp(ext, 'eps')
                ptype = '-depsc';
            elseif strcmp(ext, 'tif') || strcmp(ext, 'tiff')
                ptype = '-dtiff';
            elseif strcmp(ext, 'jpg') || strcmp(ext, 'jpeg')
                ptype = '-djpeg';
            elseif strcmp(ext, 'png')
                ptype = '-dpng';
            else
                error('File extension not recognized.')
            end
                        
            fname = [obj.conditions{i} '_' utypes{j} '.' ext];
            fname = [folder filesep strjoin(strsplit(fname, ':'), '_')];
            print(ptype, fname)
            
            close
        end
    end

end

