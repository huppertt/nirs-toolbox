function h = plotmesh( nodes, faces, values, vmax, thresh, cmap,axis_handle)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
if(nargin<7)
    axis_handle=gca;
end
    
    if nargin < 3 || isempty(values)
        h = patch(axis_handle,'vertices',nodes,'faces',faces, ...
            'facecolor',[0.9 0.9 0.9], ...
            'edgecolor','none'); 
        camlight; %(axis_handle); 
        lighting(axis_handle,'gouraud');
    else
    
        if nargin < 5 || isempty(thresh)
            thresh = 0;
        end

        if nargin < 4 || isempty(vmax)
            lst=find(abs(values)~=Inf & ~isnan(values));
            vmax = ceil( max(abs(values(lst))) );
        end

        if nargin < 6 || isempty(cmap)
            try
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',5001) )');
            catch
                cmap = parula(5001);
            end
        end

        n = size(cmap,1);
        
        cmap(round((n-1)/2+1),:) = 0.9;

        z = linspace(-vmax, vmax, n);
        lst = abs(z) < thresh;

        cmap(lst,:) = repmat( cmap(round((n-1)/2+1), :), [sum(lst) 1] );

        values(isnan(values))=0;

        h = patch(axis_handle,'vertices',nodes,'faces',faces, ...
            'FaceVertexCdata',values, ...
            'facecolor','interp', ...
            'edgecolor','none'); 
        camlight(axis_handle); 
        lighting(axis_handle,'gouraud');
        colormap(cmap), colorbar('peer',axis_handle);
        caxis(axis_handle,[-vmax vmax]);
        daspect(axis_handle,[1 1 1]);
        axis(axis_handle,'off');
      
    end


end

