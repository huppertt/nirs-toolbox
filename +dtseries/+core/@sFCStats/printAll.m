function printAll( obj, vtype, vrange, thresh, folder, ext, flip )
    %% PRINTALL - draws and saves figures to a folder
    % 
    % Args:
    %     vtype, vrange, thresh  -  same as draw function
    %     folder  -  directory to save to
    %     ext     -  file extension of images (eps, tif, or jpg)

    if(nargin<7 || isempty(flip))
    flip=[1 1];
end



if ~exist(folder, 'dir')
    mkdir(folder);
end


h=obj.draw(vtype, vrange, thresh,flip);
for i = 1:length(h)
    f=h(i);
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
    
    fname = fullfile(folder,[f.Name '.' ext]);
    print(f,ptype, fname)
    
    close(f);
    
end
