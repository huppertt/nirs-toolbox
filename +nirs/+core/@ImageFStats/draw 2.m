function draw( obj, fmax, thresh)

    %% draw - Draws channelwise F-statisitcs on a probe.
    % Args:
    %     fmax    - display range from 0 to fmax
    %     thresh  - either a scalar such that F-stats > thresh are significant or
    %               a string specifying statistical significance (e.g. 'p < 0.05')
    % 
    % Examples:
    %     stats.draw( 10, 'q < 0.1' )
    %     stats.draw( 10, 5 )
    
    % significance masking
    if nargin < 3
        mask = ones(size(obj.F)) > 0;
    
    elseif isscalar(thresh)
        mask = abs(obj.F) > thresh;
        
    elseif isvector(thresh) && isnumeric(thresh)
        mask = obj.F < thresh(1) | obj.F > thresh(2);
        
    elseif isstr(thresh)
        % takes in string in form of 'p < 0.05' or 'q < 0.10'
        s = strtrim( strsplit( thresh, '<' ) );
        
        mask = obj.(s{1}) < str2double(s{2});
    end
    
    % display range
    if nargin < 2 || isempty(fmax)
        fmax = max(obj.F);
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
    [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',4001) )');
    cmap = cmap(2000:end,:);
    
    z = linspace(0, fmax, size(cmap,1))';
    
    for iCond = 1:length( uconds )
        for iType = 1:length(utypes)
            lst = strcmp( types, utypes(iType) ) & ...
                strcmp( obj.variables.cond, uconds(iCond) );
            
            % values
            vals = obj.F(lst);
            
            % this mask
            m = mask(lst);
            
            % map to colors
            idx = bsxfun(@minus, vals', z);
            [~, idx] = min(abs(idx), [], 1);
            
            colors = cmap(idx, :);
            
            figure;
            obj.mesh.draw(vals.*m,fmax,0,cmap);
            caxis([0 fmax]);
%             % line styles
%             lineStyles = {};
%             for i = 1:length(idx)
%                 if m(i)
%                     lineStyles(i,:) = {'LineStyle', '-', 'LineWidth', 8};
%                 else
%                     lineStyles(i,:) = {'LineStyle', '--', 'LineWidth', 4};
%                 end
%             end
%             
%             figure;
%             obj.probe.draw(colors, lineStyles);
%             c = colorbar; colormap(cmap); caxis([0 fmax]);
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
end
