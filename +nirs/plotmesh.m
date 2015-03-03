function h = plotmesh( nodes, faces, values )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    
    if nargin < 3
        h = patch('vertices',nodes,'faces',faces, ...
            'facecolor',[0.8 0.8 0.8], ...
            'edgecolor','none'); 
        camlight; lighting gouraud
    else
        try
            cmap = flipud( cbrewer('div','RdBu',5001) );
        catch
            cmap = parula(5001);
        end

        cmax = max(abs(values));
        h = patch('vertices',nodes,'faces',faces, ...
            'FaceVertexCdata',values, ...
            'facecolor','interp', ...
            'edgecolor','none'); 
        camlight; lighting gouraud
        colormap(cmap), colorbar
        caxis([-cmax cmax])
    end


end

