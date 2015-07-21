function h = plotmesh( nodes, faces, values, vmax, thresh, cmap )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    
    if nargin < 3 || isempty(values)
        h = patch('vertices',nodes,'faces',faces, ...
            'facecolor',[0.9 0.9 0.9], ...
            'edgecolor','none'); 
        camlight; lighting gouraud
    else
    
        if nargin < 5 || isempty(thresh)
            thresh = 0;
        end

        if nargin < 4 || isempty(vmax)
            vmax = ceil( max(abs(values)) );
        end

        if nargin < 6 || isempty(cmap)
            try
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',5001) )');
            catch
                cmap = parula(5001);
            end
        end

        n = size(cmap,1);
        
        cmap((n-1)/2+1,:) = 0.9;

        z = linspace(-vmax, vmax, n);
        lst = abs(z) < thresh;

        cmap(lst,:) = repmat( cmap((n-1)/2+1, :), [sum(lst) 1] );


        h = patch('vertices',nodes,'faces',faces, ...
            'FaceVertexCdata',values, ...
            'facecolor','interp', ...
            'edgecolor','none'); 
        camlight; lighting gouraud
        colormap(cmap), colorbar
        caxis([-vmax vmax])
        
    end


end

