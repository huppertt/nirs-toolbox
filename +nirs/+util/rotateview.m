function rotateview(h,whichview)
% Rotates the 3D view of a figure;
% Inputs:
%   h- handles of the axis holding the mesh object


switch(whichview)
    case('left')
        viewfcn={@(h)view(h,-90,0); ...  % Change view to the left
            @(h)camroll(h,0); ... % rotate the image by 90degress
            @(h)delete(findobj('type','light','parent',h));...
            @(h)light('parent',h,'position',get(h,'cameraPosition'));...
            @(h)axis(h,'off')}; % Turn off the axis
    case('right')
        viewfcn={@(h)view(h,90,0); ...  % Change view to the right
            @(h)camroll(h,0); ... % rotate the image by 90degress
            @(h)delete(findobj('type','light','parent',h));...
            @(h)light('parent',h,'position',get(h,'cameraPosition'));...
            @(h)axis(h,'off')}; % Turn off the axis
    case{'posterior' , 'occipital' , 'back'}
        viewfcn={@(h)view(h,0,0); ...  % Change view to the left
            @(h)camroll(h,0); ... % rotate the image by 90degress
            @(h)delete(findobj('type','light','parent',h));...
            @(h)light('parent',h,'position',get(h,'cameraPosition'));...
            @(h)axis(h,'off')}; % Turn off the axis
    case{'superior' , 'top'}
        viewfcn={@(h)view(h,0,90); ...  % Change view to the top
            @(h)camroll(h,0); ... % rotate the image by 180degress
            @(h)delete(findobj('type','light','parent',h));...
            @(h)light('parent',h,'position',get(h,'cameraPosition'));...
            @(h)axis(h,'off')}; % Turn off the axis
    case{'inferior' , 'bottom'}
        viewfcn={@(h)view(h,0,-90); ...  % Change view to the left
            @(h)camroll(h,180); ... % rotate the image by 90degress
            @(h)delete(findobj('type','light','parent',h));...
            @(h)light('parent',h,'position',get(h,'cameraPosition'));...
            @(h)axis(h,'off')}; % Turn off the axis
    otherwise
        viewfcn={@(h)view(h,180,0); ...  % Change view to the left
            @(h)camroll(h,0); ... % rotate the image by 90degress
            @(h)delete(findobj('type','light','parent',h));...
            @(h)light('parent',h,'position',get(h,'cameraPosition'));...
            @(h)axis(h,'off')}; % Turn off the axis
end

cellfun(@(f)f(h),viewfcn,'UniformOutput', false);
return