function h = draw( obj, vtype, vrange, thresh ,powerthresh,viewpt)
%% draw - Draws channelwise values on a probe.
% Args:
%     vtype   - either 'beta' or 'tstat'
%     vrange  - range to display; either a scalar or vector with 2 elements
%     thresh  - either a scalar such that values > thresh are significant or
%               a string specifying statistical significance (e.g. 'p < 0.05')
%     pwerthresh- (e.g. beta>.8) applies a mask based on type-II error
%     viewpt - viewing position {['frontal' -default],'left','right','posterior','superior','inferior'}
%
% Examples:
%     stats.draw( 'tstat', [-5 5], 'q < 0.1' )
%     stats.draw( 'tstat', 5, 3 )

% type is either beta or tstat
if nargin < 2, vtype = 'tstat'; end

values = obj.(vtype);

% range to show
if nargin < 3 || isempty(vrange)
    vmax    = max(abs(values(:)));
    vrange  = vmax*[-1 1];
end

% significance mask
if nargin < 4 || isempty(thresh)
    mask = ones(size(values)) > 0;
    thresh='p<1';
    
elseif isscalar(thresh)
    mask = abs(values) > thresh;
    thresh=['p<' num2str(2*tcdf(-abs(thresh), obj.dfe))];
elseif isvector(thresh) && isnumeric(thresh)
    mask = values < thresh(1) | values > thresh(2);
    thresh=['p<' num2str(2*tcdf(-abs(thresh(2)), obj.dfe))]
elseif isstr(thresh)
    % takes in string in form of 'p < 0.05' or 'q < 0.10'
    s = strtrim( strsplit( thresh, '<' ) );
    
    mask = obj.(s{1}) < str2double(s{2});
end

if(nargin>=5) && ~isempty(powerthresh)
    pwr = obj.power(thresh);
    s = strtrim( strsplit( powerthresh, '>' ) );
    mask = mask.*(pwr > str2double(s{2}));
else
    % Leave mask alone
end
if(nargin<6)
    viewpt={'frontal'};
end


% meas types
types = obj.variables.type;

% convert to strings for consistency in loop below
if any(isnumeric(types))
    types = cellfun(@(x) {num2str(x)}, num2cell(types));
end

% unique types
utypes = unique(types, 'stable');

% unique conditions
uconds = unique(obj.variables.cond, 'stable');

% colormap
[~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
z = linspace(vrange(1), vrange(2), size(cmap,1))';


if(~iscell(viewpt))
    viewpt={viewpt};
end

for i=1:length(viewpt)
    if(isa(viewpt{i}, 'function_handle'))
        viewfcn{i}=viewpt{i};
    elseif(~isempty(viewpt{i}))
        switch(lower(viewpt{i}))
            case('left')
                viewfcn{i}={@(h)view(h,-90,0); ...  % Change view to the left
                    @(h)camroll(h,0); ... % rotate the image by 90degress
                    @(h)delete(findobj('type','light','parent',h));...
                    @(h)light('parent',h,'position',get(h,'cameraPosition'));...
                    @(h)axis(h,'off')}; % Turn off the axis
            case('right')
                viewfcn{i}={@(h)view(h,90,0); ...  % Change view to the right
                    @(h)camroll(h,0); ... % rotate the image by 90degress
                    @(h)delete(findobj('type','light','parent',h));...
                    @(h)light('parent',h,'position',get(h,'cameraPosition'));...
                    @(h)axis(h,'off')}; % Turn off the axis
            case{'posterior' , 'occipital' , 'back'}
                viewfcn{i}={@(h)view(h,0,0); ...  % Change view to the left
                    @(h)camroll(h,0); ... % rotate the image by 90degress
                    @(h)delete(findobj('type','light','parent',h));...
                    @(h)light('parent',h,'position',get(h,'cameraPosition'));...
                    @(h)axis(h,'off')}; % Turn off the axis
            case{'superior' , 'top'}
                viewfcn{i}={@(h)view(h,0,90); ...  % Change view to the top
                    @(h)camroll(h,0); ... % rotate the image by 180degress
                    @(h)delete(findobj('type','light','parent',h));...
                    @(h)light('parent',h,'position',get(h,'cameraPosition'));...
                    @(h)axis(h,'off')}; % Turn off the axis
            case{'inferior' , 'bottom'}
                viewfcn{i}={@(h)view(h,0,-90); ...  % Change view to the left
                    @(h)camroll(h,180); ... % rotate the image by 90degress
                    @(h)delete(findobj('type','light','parent',h));...
                    @(h)light('parent',h,'position',get(h,'cameraPosition'));...
                    @(h)axis(h,'off')}; % Turn off the axis
            otherwise
                viewfcn{i}={@(h)view(h,180,0); ...  % Change view to the left
                    @(h)camroll(h,0); ... % rotate the image by 90degress
                    @(h)delete(findobj('type','light','parent',h));...
                    @(h)light('parent',h,'position',get(h,'cameraPosition'));...
                    @(h)axis(h,'off')}; % Turn off the axis
        end
    else
         viewfcn{i}={@(h)view(h,180,0); ...  % Change view to the left
                    @(h)camroll(h,0); ... % rotate the image by 90degress
                    @(h)delete(findobj('type','light','parent',h));...
                    @(h)light('parent',h,'position',get(h,'cameraPosition'));...
                    @(h)axis(h,'off')}; % Turn off the axis
        
    end
end

h=[];

        

for iCond = 1:length( uconds )
    for iType = 1:length(utypes)
        lst = strcmp( types, utypes(iType) ) & ...
            strcmp( obj.variables.cond, uconds(iCond) );
        m = mask(lst);
        % values
        vals = values(lst);
        
       
        
        % map to colors
        idx = bsxfun(@minus, vals', z);
        [~, idx] = min(abs(idx), [], 1);
        
        colors = cmap(idx, :);
        
        if(isstr(thresh));
            thresh=obj.getCritT(thresh);
        end
        
        figure;
        h(end+1)=axes;
        obj.mesh.draw(vals.*m,vrange(2),thresh,cmap);
        %c = colorbar; colormap(cmap); caxis(vrange);
        %             a = gca;
        %
        %             ap = get(a, 'Position');
        %
        %             cp = get(c, 'Position');
        %             cp(3) = 0.5*cp(3);
        %             set(c, 'Position', cp);
        %             set(a, 'Position', ap);
        
        title([utypes{iType} ' : ' uconds{iCond}], 'Interpreter','none')
        
        
        
        
    end
end

ff=figure('visible','off');
for j=1:length(viewfcn)
    sp(j)=subplot(1,length(viewfcn),j);
    colorbar('SouthOutside');
end

for i=1:size(h,2)
    f=figure(get(h(1,i),'parent'));
    cb=colorbar;
    
    set(h(1,i),'position',get(sp(1),'Position'));
    set(cb,'AxisLocation','in','Location','southoutside');
    for j=2:length(viewfcn)
        h(j,i)=axes('Position',get(sp(j),'Position'),'parent',f);
        copyobj(get(h(1,i),'children'),h(j,i));
        set(h(j,i),'clim',get(h(1,i),'clim'));
        axes(h(j,i));
        cb=colorbar('AxisLocation','in','Location','southoutside');
        set(h(j,i),'position',get(sp(j),'Position'));
        title(h(j,i),get(get(h(1,i),'title'),'string'),'Interpreter','none')
    end
    axis(h(:,i),'off');
    set(f,'color','w');
    
end

close(ff);

for idx=1:length(viewfcn)
    cellfun(@(f)arrayfun(f,handle(h(idx,:)),'UniformOutput', false),viewfcn{idx},'UniformOutput', false);
end
end