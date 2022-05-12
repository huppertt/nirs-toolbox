function f = draw( obj, vtype, vrange, thresh,figH,ConditionsShown,TypesShown)
%% draw - Draws channelwise values on a probe.
% Args:
%     vtype   - either 'beta' or 'tstat'
%     vrange  - range to display; either a scalar or vector with 2 elements
%     thresh  - either a scalar such that values > thresh are significant or
%               a string specifying statistical significance (e.g. 'p < 0.05')
%
% Examples:
%     stats.draw( 'tstat', [-5 5], 'q < 0.1' )
%     stats.draw( 'tstat', 5, 3 )

% type is either beta or tstat
if nargin < 2, vtype = 'tstat'; end

% Need to make sure variables and link tables are in same order
obj = obj.sorted();

v=[];
for ii=1:length(obj)
    values{ii} = obj(ii).(vtype);
    v=[v; values{ii}(:)];
end


% range to show
if nargin < 3 || isempty(vrange)
    vmax    = max(abs(v(:)));
    vrange  = vmax*[-1 1];
end



for ii=1:length(values)
    % significance mask
    if nargin < 4
        mask{ii} = ones(size(values{ii})) > 0;
        
    elseif isscalar(thresh)
        
        mask{ii} = abs(values{ii}) > thresh;
        
    elseif isvector(thresh) && isnumeric(thresh)
        mask{ii} = values{ii} < thresh(1) | values{ii} > thresh(2);
        
    elseif isstr(thresh)
        % takes in string in form of 'p < 0.05' or 'q < 0.10'
        s = strtrim( strsplit( thresh, '<' ) );
        
        mask{ii} = obj(ii).(s{1}) < str2double(s{2});
    end
end

visible='on';
callers = dbstack(1);
if ~isempty(callers)
    if strcmp(callers(1).name,'printAll')
        visible='off';
    end
end

% meas types
if(length(obj)>1)
    types=vertcat(obj.variables);
    types=types.type;
else
    types = obj.variables.type;
end
% convert to strings for consistency in loop below
if any(isnumeric(types))
    types = cellfun(@(x) {num2str(x)}, num2cell(types));
end


% unique types
utypes = unique(types, 'stable');
% unique conditions


if(nargin>6 && ~isempty(TypesShown))
    if(ischar(TypesShown))
        TypesShown={TypesShown};
    end
    utypes={utypes{ismember(utypes,TypesShown)}};
end

if(length(obj)>1)
    cond=vertcat(obj.variables);
    cond=cond.cond;
    uconds = unique(cond, 'stable');
else
    uconds = unique(obj.variables.cond, 'stable');
end

if(nargin>5 && ~isempty(ConditionsShown))
     if(ischar(ConditionsShown))
        ConditionsShown={ConditionsShown};
    end
     uconds={ uconds{ismember( uconds,ConditionsShown)}};
end


% colormap
[~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',128) )');
z = linspace(vrange(1), vrange(2), size(cmap,1))';

hind = 0;
for iType = 1:length(utypes)
    
    if(isa(obj(1).probe,'nirs.core.ProbeROI'))
        hind = hind + 1;
        
        if(nargin<5 || isempty(figH))
            f(hind)=figure('Visible',visible);
        else
            f(hind)=figH;
        end
        
        set(f(hind),'name',[utypes{iType}]);
        for ii=1:length(obj)
            
            
            if(length(obj)>1)
                n=ceil(sqrt(length(obj)));
                a=subplot(n,n,ii);
            else
                a = get(f(hind),'CurrentAxes');
            end
            types = obj(ii).variables.type;
            if any(isnumeric(types))
                types = cellfun(@(x) {num2str(x)}, num2cell(types));
            end
            
            lst = strcmp( types, utypes(iType));
            if(isempty(lst))
                continue;
            end
            % values
            vals = values{ii}(lst);
            
            % this mask
            m = mask{ii}(lst);
            
            % map to colors
            idx = bsxfun(@minus, vals', z);
            [~, idx] = min(abs(idx), [], 1);
            
            colors = cmap(idx, :);
            
            % line styles
            lineStyles = {};
            for i = 1:length(idx)
                if m(i)
                    lineStyles(i,:) = {'LineStyle', '-', 'LineWidth', 8};
                else
                    lineStyles(i,:) = {'LineStyle', '--', 'LineWidth', 4};
                end
            end
            
            
            if(isempty(a))
                figure(f(hind))
                a=axes;
            end
            
           
            obj(ii).probe.draw(colors, lineStyles,a,obj(ii).variables.cond(lst),vals);
            
            if(~strcmp(lower(obj(ii).probe.defaultdrawfcn),'bar') & ...
                    ~(isempty(obj(ii).probe.defaultdrawfcn) & ...
                    isa(obj(ii).probe,'nirs.core.ProbeROI')))
                c = colorbar; colormap(cmap); caxis(vrange);
                ap = get(a, 'Position');
                
                cp = get(c, 'Position');
                cp(3) = 0.5*cp(3);
                set(c, 'Position', cp);
                set(a, 'Position', ap);
            else
                set(a,'Ylim',vrange);
            end
            
            title(utypes{iType}, 'Interpreter','none')
            
        end
        
    else
        
        for iCond = 1:length( uconds )
            
            hind = hind + 1;
            
            if(nargin<5 || isempty(figH))
                f(hind)=figure('Visible',visible);
            else
                f(hind)=figH;
            end
            
            set(f(hind),'name',[utypes{iType} ' : ' uconds{iCond}]);
            
            for ii=1:length(obj)
                
                if(length(obj)>1)
                    n=ceil(sqrt(length(obj)));
                    a=subplot(n,n,ii);
                else
                    a = get(f(hind),'CurrentAxes');
                end
                
                types = obj(ii).variables.type;
                if any(isnumeric(types))
                    types = cellfun(@(x) {num2str(x)}, num2cell(types));
                end
                
                lst = strcmp( types, utypes(iType) ) & ...
                    strcmp( obj(ii).variables.cond, uconds(iCond) );
                if(isempty(lst))
                    continue;
                end
                % values
                vals = values{ii}(lst);
                
                % this mask
                m = mask{ii}(lst);
                
                % map to colors
                idx = bsxfun(@minus, vals', z);
                [~, idx] = min(abs(idx), [], 1);
                
                colors = cmap(idx, :);
                
                % line styles
                lineStyles = {};
                for i = 1:length(idx)
                    if m(i)
                        lineStyles(i,:) = {'LineStyle', '-', 'LineWidth', 8};
                    else
                        lineStyles(i,:) = {'LineStyle', '--', 'LineWidth', 4};
                    end
                end
                
                
                if(isempty(a))
                    figure(f(hind))
                    a=axes;
                end
                hh=obj(ii).probe.draw(colors, lineStyles,a);
                
                set(hh(isnan(vals)),'visible','off');
                
                c = colorbar; colormap(cmap); caxis(vrange);
                
                if(all(get(gcf,'color')==0))
                   set(c,'Color','w')
                end
                ap = get(a, 'Position');
                
                cp = get(c, 'Position');
                cp(3) = 0.5*cp(3);
                set(c, 'Position', cp);
                set(a, 'Position', ap);
                
                title([utypes{iType} ' : ' uconds{iCond}], 'Interpreter','none')
            end
        end
    end
end

end