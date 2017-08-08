classdef TreeDrawer
%TreeDrawer Class to draw decision tree. The code has been copied from
%classregtree/view with minor changes.

%   Copyright 2012-2014 The MathWorks, Inc.
    
    methods(Static)
        function fig = setupfigure(doclass)
            %SETUPFIGURE Set up uicontrols on decision tree figure.
            
            fig = figure('IntegerHandle','off', 'NumberTitle','off', ...
                'Units','points','PaperPositionMode','auto',...
                'Tag','tree viewer');
            ax = axes('Parent',fig,'UserData',cell(1,4),'XLim',0:1,'YLim',0:1);
            
            % Set default print options
            pt = printtemplate;
            pt.PrintUI = 0;
            set(fig,'PrintTemplate',pt)
            
            if doclass
                figtitle = getString(message('stats:classregtree:view:ClassificationTreeViewer'));
            else
                figtitle = getString(message('stats:classregtree:view:RegressionTreeViewer'));
            end
            
            % Arrange figure contents
            pos = [0 0 1 1];
            set(ax,'Visible','off','XLim',0:1,'YLim',0:1,'Position',pos);
            set(ax,'Units','points');
            apos = get(ax,'Position');
            fpos = get(fig,'Position');
            hframe = uicontrol(fig,'Units','points','Style','frame',...
                'Position',[0 0 1 1],'Tag','frame');
            
            % Tip-related items
            h=uicontrol(fig,'units','points','Tag','clicktext',...
                'String',getString(message('stats:classregtree:view:ClickToDisplay')),...
                'style','text','HorizontalAlignment','left','FontWeight','bold');
            extent = get(h,'Extent');
            theight = extent(4);
            aheight = apos(4);
            tbottom = aheight - 1.5*theight;
            posn = [2, tbottom, 150, theight];
            set(h,'Position',posn);
            textpos = posn;
            e = get(h,'Extent');
            if doclass
                choices = ...
                    { getString(message('stats:classregtree:view:NodeIdentity')), ...
                    getString(message('stats:classregtree:view:NodeVariableRanges')), ...
                    getString(message('stats:classregtree:view:NodeClassMembership')), ...
                    getString(message('stats:classregtree:view:NodeEstimatedProbabilities')) };
            else
                choices = ...
                    { getString(message('stats:classregtree:view:NodeIdentity')),...
                    getString(message('stats:classregtree:view:NodeVariableRanges')),...
                    getString(message('stats:classregtree:view:NodeStatistics')) };
            end
            strlengths = cellfun('length',choices);
            [~,longest] = max(strlengths);
            h=uicontrol(fig,'units','points','position',[0 0 1 1],'Tag','clicklist',...
                'String',choices{longest}, 'Style','pop','BackgroundColor',ones(1,3),...
                'Callback',@classreg.learning.treeutils.TreeDrawer.removelabels);            
            hext = get(h,'Extent');
            posn = [e(1)+e(3)+2, aheight-1.25*theight, hext(3)+40, theight];
            set(h,'Position',posn);
            set(h,'String',choices);
            set(ax,'Position',[0 0 apos(3) tbottom]);
            set(fig,'Toolbar','figure', 'Name',figtitle,'HandleVisibility','callback');
            
            % Magnification items
            textpos(1) = posn(1) + posn(3) + 10;
            h=uicontrol(fig,'units','points','Tag','magtext','Position',textpos,...
                'String',getString(message('stats:classregtree:view:PlotMagnification')),...
                'style','text','HorizontalAlignment','left','FontWeight','bold');
            e = get(h,'Extent');
            textpos(3) = e(3);
            set(h,'Position',textpos);
            h=uicontrol(fig,'units','points','position',[0 0 1 1],'Tag','maglist',...
                'String','x', 'Style','pop','BackgroundColor',ones(1,3),...
                'Callback',@domagnif);
            adjustcustomzoom(h,false);
            hext = get(h,'Extent');
            posn = [textpos(1)+textpos(3)+2, posn(2), hext(3)+80, posn(4)];
            set(h,'Position',posn);
            
            % Prune-related items
            textpos(1) = posn(1) + posn(3) + 10;
            h = uicontrol(fig,'units','points','position',textpos,'Tag','prunelabel',...
                'Style','text','HorizontalAlignment','left',...
                'FontWeight','bold',...
                'String',getString(message('stats:classregtree:view:PruneLevel')));
            e = get(h,'Extent');
            textpos(3) = e(3);
            set(h,'Position',textpos);
            
            posn(1) = textpos(1) + textpos(3) + 5;
            posn(2) = posn(2) - 0.25*e(4);
            posn(3) = 60;
            posn(4) = 1.5*e(4);
            textpos(1) = posn(1) + 3;
            textpos(3) = posn(3) - 6;
            uicontrol(fig,'Style','frame','Units','points','Position',posn,...
                'Tag','pruneframe');
            uicontrol(fig,'units','points','position',textpos,'Tag','prunelev',...
                'Style','text','HorizontalAlignment','left',...
                'FontWeight','bold','String','1234 of 9999');
            
            % Create an arrow for labeling button controls
            fcolor = get(fig,'Color');
            ar = ...
                [1 1 1 1 1 1 1 1 1 1 1 1 1
                1 1 1 1 1 1 1 1 1 1 1 1 1
                1 0 0 0 0 0 0 0 0 0 0 0 1
                1 1 0 0 0 0 0 0 0 0 0 1 1
                1 1 1 0 0 0 0 0 0 0 1 1 1
                1 1 1 1 0 0 0 0 0 1 1 1 1
                1 1 1 1 1 0 0 0 1 1 1 1 1
                1 1 1 1 1 1 0 1 1 1 1 1 1
                1 1 1 1 1 1 1 1 1 1 1 1 1
                1 1 1 1 1 1 1 1 1 1 1 1 1];
            ar = repmat(ar,[1 1 3]);
            ar(:,:,1) = min(ar(:,:,1),fcolor(1));
            ar(:,:,2) = min(ar(:,:,2),fcolor(2));
            ar(:,:,3) = min(ar(:,:,3),fcolor(3));
            
            posn(1) = posn(1) + posn(3);
            posn(4) = posn(4)/2;
            posn(3) = posn(4);
            uicontrol(fig,'units','points','position',posn,'Tag','prune',...
                'CData',ar(end:-1:1,:,:),...
                'Style','pushbutton','Callback',@growprune);
            posn(2) = posn(2) + posn(4);
            uicontrol(fig,'units','points','position',posn,'Tag','grow',...
                'CData',ar,...
                'Style','pushbutton','Callback',@growprune);
            
            if posn(1)+posn(3) > fpos(3)
                fpos(3) = posn(1)+posn(3);
                set(fig,'Position',fpos);
                apos(3) = fpos(3);
                set(ax,'Position',apos);
            end
            
            % Adjust frame position
            lowest = min(posn(2),textpos(2))-2;
            frpos = apos;
            frpos(4) = 1.1*(apos(4)-lowest);
            frpos(2) = apos(4) - frpos(4);
            set(hframe,'Position',frpos);
            
            % Add scroll bars, start out invisible
            h1 = uicontrol(fig,'Style','slider','Tag','hslider','Visible','off',...
                'Units','points','Callback',@dopan);
            p1 = get(h1,'Position');
            sw = p1(4);               % default slider width
            p1(1:2) = 1;
            p1(3) = fpos(3)-sw;
            set(h1,'Position',p1);
            p2 = [fpos(3)-sw, sw, sw, frpos(2)-sw];
            uicontrol(fig,'Style','slider','Tag','vslider','Visible','off',...
                'Units','points','Position',p2,'Callback',@dopan);
            
            % Add new menu before the Window menu
            hw = findall(fig,'Type','uimenu','Tag','figMenuWindow');
            h0 = uimenu(fig,'Label',getString(message('stats:classregtree:view:TreeMenu')),...
                'Position',get(hw,'Position'));
            uimenu(h0, 'Label',getString(message('stats:classregtree:view:TreeMenu_ShowFullTree')), ...
                'Position',1,'Tag','menufull','Checked','on','Callback',@domenu);
            uimenu(h0, 'Label',getString(message('stats:classregtree:view:TreeMenu_ShowUnprunedNodes')), ...
                'Position',2,'Tag','menuunpr','Checked','off','Callback',@domenu);
            uimenu(h0, 'Label',getString(message('stats:classregtree:view:TreeMenu_LabelBranchNodes')), ...
                'Position',3,'Tag','menubr','Checked','on','Callback',@domenu,'Separator','on');
            uimenu(h0, 'Label',getString(message('stats:classregtree:view:TreeMenu_LabelLeafNodes')), ...
                'Position',4,'Tag','menuleaf','Checked','on','Callback',@domenu);
            
            set(fig,'ResizeFcn',@resize);
        end
        
        function adjustmenu(fig,helpviewer)
            %ADJUSTMENU Adjust contents of curve fit plot menus and toolbar
            
            % Remove some menus entirely
            badTags = {'figMenuEdit' 'figMenuView' 'figMenuInsert'};
            h = findall(fig, 'Type','uimenu', 'Parent',fig);
            tagFun = @(x) get(x,'Tag');
            foundTags = arrayfun(tagFun,h,'UniformOutput',false);
            tf = ismember(foundTags,badTags);
            delete(h(tf));
            h(tf) = [];
            
            % Add or remove some items from other menus
            % Fix FILE menu
            h0 = findall(h, 'Type','uimenu', 'Tag','figMenuFile');
            h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
            badTags = {'printMenu' 'figMenuFilePrintPreview'};
            foundTags = arrayfun(tagFun,h1,'UniformOutput',false);
            tf = ismember(foundTags,badTags);
            delete(h1(tf));
            
            % Fix TOOLS menu
            h0 = findall(h, 'Type','uimenu', 'Tag','figMenuTools');
            h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
            badTags = {'figMenuZoomIn' 'figMenuZoomOut' 'figMenuPan'};
            foundTags = arrayfun(tagFun,h1,'UniformOutput',false);
            tf = ismember(foundTags,badTags);
            delete(h1(~tf));
            
            % Fix HELP menu
            h0 = findall(h, 'Type','uimenu', 'Tag','figMenuHelp');
            h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
            delete(h1);
            uimenu(h0, 'Label', getString(message('stats:classregtree:view:HelpTreeViewer')), ...
                'Position',1,'Callback',helpviewer);
           
            % Remove icons that don't apply here.  Keep zoom, PAN and print only.
            h0 = findall(fig,'Type','uitoolbar');
            h1 = findall(h0,'Parent',h0);
            badTags = {'Exploration.Pan' 'Exploration.ZoomOut' 'Exploration.ZoomIn'};
            foundTags = arrayfun(tagFun,h1,'UniformOutput',false);
            tf = ismember(foundTags,badTags);
            delete(h1(~tf));
        end

        function [X,Y] = drawtree(tree,fig,nodevalue,varnames,curlevel,classnames)
            %DRAWTREE Draw decision tree into specified figure.
            
            ax = get(fig,'CurrentAxes');
            splitvar = tree.CutVar;
            cutpoint = tree.CutPoint;
            cutcateg = tree.CutCategories;
            parentnode = tree.Parent;
            nonroot = parentnode~=0;
            
            % Get coordinates of nodes
            isleaf = splitvar==0;
            [X,Y] = layouttree(tree,isleaf);
            
            % Which leaves and branches are pruned?
            if ~isempty(tree.PruneList)
                if curlevel==0
                    isbranch = tree.IsBranch;
                else
                    isbranch = (tree.PruneList>curlevel);
                end
            else
                isbranch = ~isleaf;
            end
            if any(isbranch)
                isleaf = false(size(isbranch));
                c = tree.Children(isbranch,:);
                c = c(c>0);
                isleaf(c) = 1;
                isleaf = isleaf & ~isbranch;
            else
                isleaf = ~nonroot;
            end
            pruned = ~(isleaf | isbranch);
            
            branchnodes = find(isbranch);
            leafnodes = find(isleaf);
            
            % Get coordinates of connecting lines
            p = parentnode(nonroot & ~pruned);
            x = [X(nonroot & ~pruned)'; X(p)'; NaN+p'];
            y = [Y(nonroot & ~pruned)'; Y(p)'; NaN+p'];
            
            % Plot nodes and connections for unpruned nodes, but stop listening to axis
            % and remember some things that may get changed during the plotting
            axislistener(ax,false);
            xlim = get(ax,'XLim');
            ylim = get(ax,'YLim');
            ud = get(ax,'UserData');
            h = plot(X(branchnodes),Y(branchnodes),'b^', ...
                X(leafnodes),Y(leafnodes),'b.', ...
                x(:),y(:),'b-','Parent',ax);
            set(ax,'UserData',ud,'Visible','on','XLim',xlim,'YLim',ylim);
            set(ax,'Color','none');
            set(ax,'XTick',[]);
            set(ax,'YTick',[]);
            axislistener(ax,true);
            
            % Same for pruned nodes
            t = nonroot & pruned;
            p = parentnode(t);
            x = [X(t)'; X(p)'; NaN+p'];
            y = [Y(t)'; Y(p)'; NaN+p'];
            line('Parent',ax,'XData',X(pruned),'YData',Y(pruned),'Tag','prunednode',...
                'Marker','o','Color',[.2 .2 .2],'Linestyle','none','HitTest','off',...
                'PickableParts','none');
            line('Parent',ax,'XData',x(:),'YData',y(:),'Tag','prunedconnection',...
                'Marker','none','LineStyle',':','Color',[.2 .2 .2],'HitTest','off',...
                'PickableParts','none');
            if length(h)==3
                set(h(1),'ButtonDownFcn',@labelpoint,'Tag','branch','MarkerSize',10);
                set(h(2),'ButtonDownFcn',@labelpoint,'Tag','leaf','MarkerSize',20);
                set(h(end),'HitTest','off','PickableParts','none','Tag','connection');
            else
                set(h,'ButtonDownFcn',@labelpoint,'Tag','leaf','MarkerSize',20);
            end
            
            % Label leaf nodes with class, branch nodes with split rule
            if ~isempty(classnames)
                ctext = classnames(nodevalue(leafnodes));
            else
                ctext = num2str(nodevalue(leafnodes));
            end
            
            h = findobj(fig,'Tag','menuleaf');
            vis = get(h,'Checked');
            text(X(leafnodes),Y(leafnodes),ctext,'HitTest','off','PickableParts','none','Parent',ax,...
                'VerticalAlignment','top','HorizontalAlignment','center',...
                'Tag','leaflabel','Clipping','on','Visible',vis,'Interpreter','none');
            
            lefttext = cell(length(branchnodes),1);
            righttext = cell(length(branchnodes),1);
            for j=1:length(branchnodes)
                k = branchnodes(j);
                %cut = cutoff{k};
                varname = varnames{splitvar(k)};
                if ~isnan(cutpoint(k)) % numeric split
                    lefttext{j} = sprintf('%s < %g   ',varname,cutpoint(k));
                    righttext{j} = sprintf('  %s >= %g',varname,cutpoint(k));
                else
                    cats = cutcateg{k,1};
                    if length(cats)==1
                        lefttext{j} = sprintf('%s = %s   ',varname,num2str(cats,'%g '));
                    else
                        lefttext{j} = sprintf('%s %s (%s)   ',varname,getString(message('stats:classregtree:disp:ElementInSet')),num2str(cats,'%g '));
                    end
                    cats = cutcateg{k,2};
                    if length(cats)==1
                        righttext{j} = sprintf('   %s = %s',varname,num2str(cats,'%g '));
                    else
                        righttext{j} = sprintf('   %s %s (%s)',varname,getString(message('stats:classregtree:disp:ElementInSet')),num2str(cats,'%g '));
                    end
                end
            end
            
            h = findobj(fig,'Tag','menubr');
            vis = get(h,'Checked');
            text(X(branchnodes),Y(branchnodes),lefttext,'HitTest','off','PickableParts','none','Parent',ax,...
                'Tag','branchlabel','Clipping','on','Visible',vis,'Interpreter','none',...
                'HorizontalAlignment','right');
            text(X(branchnodes),Y(branchnodes),righttext,'HitTest','off','PickableParts','none','Parent',ax,...
                'Tag','branchlabel','Clipping','on','Visible',vis,'Interpreter','none',...
                'HorizontalAlignment','left');
            
            % Show pruned nodes or not as desired
            doprunegraph(fig);
            
            % Adjust axes contents
            dozoom(fig);
            
            % Adjust layout of controls to fit figure
            layoutfig(fig);
        end
        
        function updatelevel(fig,curlevel,tree)
            %UPDATELEVEL Update text display of current pruning level.
            
            if ~isempty(tree.PruneList)
                maxlevel = max(tree.PruneList);
            else
                maxlevel = 0;
            end
            h = findobj(fig,'Tag','prunelev');
            set(h,'String',sprintf('%s',getString(...
                message('stats:classregtree:view:PruneLevelOneOfMax',curlevel,maxlevel))));
            e = get(h,'Extent');
            p = get(h,'Position');
            p(3) = e(3);
            set(h,'Position',p);
        end
        
        function updateenable(fig)
            %UPDATEENABLE Update enabled/disabled status of buttons in figure
            
            ud = get(fig,'UserData');
            curlevel = ud{7};
            fulltree = ud{6};
            enableg = 'on';     % enable grow button?
            enablep = 'on';     % enable prune button?
            if isempty(fulltree.PruneList)
                enableg = 'off';
                enablep = 'off';
            else
                maxlevel = max(fulltree.PruneList);
                if curlevel >= maxlevel
                    enablep = 'off';
                end
                if curlevel <= 0;
                    enableg = 'off';
                end
            end
            set(findobj(fig,'tag','prune'),'Enable',enablep);
            set(findobj(fig,'tag','grow'),'Enable',enableg);
        end
        
        function removelabels(varargin)
            %REMOVELABELS Remove any labels remaining on the graph.
            
            f = gcbf;
            delete(findall(f,'Tag','LinescanMarker'));
            delete(findall(f,'Tag','LinescanText'));
        end
    end
    
end

function [X,Y] = layouttree(tree,isleaf)
%LAYOUTTREE Select x,y coordinates of tree elements.
n = size(tree.Children,1);
X = zeros(n,1);
Y = X;
layoutstyle = 1;

% Plot top node on one level, its children at next level, etc.
for j=1:n
    p = tree.Parent(j);
    if p>0
        Y(j) = Y(p)+1;
    end
end
if layoutstyle==1
    % Layout style 1
    % Place terminal nodes, then put parents above them
    
    % First get preliminary placement, used to position leaf nodes
    for j=1:n
        p = tree.Parent(j);
        if p==0
            X(j) = 0.5;
        else
            dx = 2^-(Y(j)+1);
            if j==tree.Children(p,1)
                X(j) = X(p) - dx;
            else
                X(j) = X(p) + dx;
            end
        end
    end
    
    % Now make leaf nodes equally spaced, preserving their order
    leaves = find(isleaf);
    nleaves = length(leaves);
    [~,b] = sort(X(leaves));
    X(leaves(b)) = (1:nleaves) / (nleaves+1);
    
    % Position parent nodes above and between their children
    for j=max(Y):-1:0
        a = find(~isleaf & Y==j);
        c = tree.Children(a,:);
        X(a) = (X(c(:,1))+X(c(:,2)))/2;
    end
else
    % Layout style 2
    % Spread out the branch nodes, somewhat under their parent nodes
    X(Y==0) = 0.5;
    for j=1:max(Y)
        vis = (Y==j);                  % real nodes at this level
        invis = (Y==(j-1) & isleaf);   % invisible nodes to improve layout
        nvis = sum(vis);
        nboth = nvis + sum(invis);
        x = [X(tree.Parent(vis))+1e-10*(1:nvis)'; X(invis)];
        [xx,xidx] = sort(x);
        xx(xidx) = 1:nboth;
        X(vis) = (xx(1:nvis) / (nboth+1));
    end
end

k = max(Y);
Y = 1 - (Y+0.5)/(k+1);
end

% ----------------------------------------------
function growprune(varargin)
%GROWPRUNE Expand or contract tree using optimal pruning sequence.

% Fetch stored information
h = gcbo;
fig = gcbf;
ud = get(fig,'UserData');
varnames = ud{4};
nodevalue = ud{5};
fulltree = ud{6};
curlevel = ud{7};
cnames = ud{8};

% Get optimal pruning sequence and current pruning level
prunelist = fulltree.PruneList;
if isequal(get(h,'Tag'),'prune')
    curlevel = min(max(prunelist),curlevel+1);
else
    curlevel = max(0,curlevel-1);
end

% Clear axes, then draw at new pruning level
ax = get(fig,'CurrentAxes');
delete(get(ax,'Children'));
[X,Y] = classreg.learning.treeutils.TreeDrawer.drawtree(fulltree,fig,nodevalue,varnames,curlevel,cnames);

% Remember everything
set(fig,'ButtonDownFcn',@classreg.learning.treeutils.TreeDrawer.removelabels, ...
    'UserData',{X Y 0 varnames nodevalue fulltree curlevel cnames});

classreg.learning.treeutils.TreeDrawer.updateenable(fig);
classreg.learning.treeutils.TreeDrawer.updatelevel(fig,curlevel,fulltree);
end

% ----------------------------------------------
function labelpoint(varargin)
%LABELPOINT Label point on tree in response to mouse click.

h = gcbo;
f = gcbf;
stype = get(f,'SelectionType');
if ~isequal(stype,'alt') && ~isequal(stype,'extend')
    classreg.learning.treeutils.TreeDrawer.removelabels;
end
t = get(h,'Tag');
if isequal(t,'branch') || isequal(t,'leaf')
    ud = get(f,'UserData');
    X             = ud{1}; % x coordinates
    Y             = ud{2}; % y coordinates
    %W            = ud{3}; % weights
    varnames      = ud{4}; % variable names
    nodevalue     = ud{5}; % node predictions
    tree          = ud{6}; % complete tree
    %curlevel     = ud{7}; % current pruning level
    cnames        = ud{8}; % class names
    
    doclass = ~isempty(cnames);
    
    splitvar = tree.CutVar;
    cutpoint = tree.CutPoint;
    cutcateg = tree.CutCategories;

    % Find closest point
    ax = get(f,'CurrentAxes');
    cp = get(ax,'CurrentPoint');
    [~,node] = min(abs(X-cp(1,1)) + abs(Y-cp(1,2)));

    uih = findobj(f,'Tag','clicklist');
    labeltype = get(uih,'Value');

    if isequal(labeltype,4) && doclass
        % Show fitted class probabilities
        P = tree.ClassProb;
        txt = getString(message('stats:classregtree:view:ClassProbabilities'));
        for j=1:size(P,2);
            txt = sprintf('%s\n%s = %.3g',txt,cnames{j},P(node,j));
        end

    elseif isequal(labeltype,3) && ~doclass
        % Show node statistics
        xbar = nodevalue;
        Nk = tree.NodeSize(node);
        if Nk > 1
            s = sqrt((tree.NodeRisk(node)./tree.NodeProb(node) * Nk) / (Nk - 1));
            txt = sprintf('N = %d\n%s = %g\n%s = %g',...
                Nk,...
                getString(message('stats:classregtree:view:TreeNodeMean')),...
                xbar(node),...
                getString(message('stats:classregtree:view:TreeNodeStandardDeviation')),...
                s);
        else
            txt = sprintf('N = %d\n%s = %g',...
                Nk,...
                getString(message('stats:classregtree:view:TreeNodeMean')),...
                xbar(node));
        end

    elseif isequal(labeltype,3) && doclass
        % Show class membership in data
        C = tree.ClassCount;
        N = tree.NodeSize(node);
        txt = sprintf('%s = %d',getString(message('stats:classregtree:view:TotalDataPoints')),N);
        for j=1:size(C,2);
            txt = sprintf('%s\n%d %s',txt,C(node,j),cnames{j});
        end

    elseif isequal(labeltype,1)
        % Get a display of the split rule at branch nodes,
        % or the majority class at leaf nodes
        if ~isequal(t,'branch')
            if doclass
                txt = sprintf('%s %d (%s)\n%s: %s',...
                    getString(message('stats:classregtree:view:TreeNode')),...
                    node,...
                    getString(message('stats:classregtree:view:TreeLeaf')),...
                    getString(message('stats:classregtree:view:PredictedClass')),...
                    cnames{nodevalue(node)});
            else
                txt = sprintf('%s %d (%s)\n%s: %g',...
                    getString(message('stats:classregtree:view:TreeNode')),...
                    node,...
                    getString(message('stats:classregtree:view:TreeLeaf')),...
                    getString(message('stats:classregtree:view:RegressionPrediction')),...
                    nodevalue(node));
            end
        elseif ~isnan(cutpoint(node)) %splitvar(node)>0
            txt = sprintf('%s %d (%s)\n%s:  %s < %g',...
                getString(message('stats:classregtree:view:TreeNode')),...
                node,...
                getString(message('stats:classregtree:view:TreeBranch')),...
                getString(message('stats:classregtree:view:SplittingRule')),...
                varnames{splitvar(node)},...
                cutpoint(node));
        else
            cut = cutcateg(node,:);
            cats = cut{1};
            if length(cats)==1
                txt = sprintf('%s %d (%s)\n%s:  %s = %s',...
                    getString(message('stats:classregtree:view:TreeNode')),...
                    node,...
                    getString(message('stats:classregtree:view:TreeBranch')),...
                    getString(message('stats:classregtree:view:SplittingRule')),...
                    varnames{splitvar(node)},...
                    num2str(cats,'%g '));
            else
                txt = sprintf('%s %d (%s)\n%s:  %s %s (%s)',...
                    getString(message('stats:classregtree:view:TreeNode')),...
                    node,...
                    getString(message('stats:classregtree:view:TreeBranch')),...
                    getString(message('stats:classregtree:view:SplittingRule')),...
                    varnames{splitvar(node)},...
                    getString(message('stats:classregtree:disp:ElementInSet')),...
                    num2str(cats,'%g '));
            end
        end
    elseif isequal(labeltype,2)
        % Get the conditions satisfied by points in this node
        if node==1
            txt = getString(message('stats:classregtree:view:RootOfTree'));
        else
            % Find limits for each variable along the path to this node
            nvars = max(splitvar(:));
            lims = cell(nvars,3);
            lims(:,1) = {-Inf};               % lower limit
            lims(:,2) = {Inf};                % upper limit
            lims(:,3) = num2cell((1:nvars)'); % variable number
            p = tree.Parent(node);
            c = node;
            while(p>0)
                leftright = 1 + (tree.Children(p,2)==c);
                vnum = splitvar(p);
                if ~isnan(cutpoint(p))
                    if isinf(lims{vnum,3-leftright})
                        lims{vnum,3-leftright} = cutpoint(p);
                    end
                else
                    if ~iscell(lims{vnum,1})
                        vcut = cutcateg(p,:);
                        lims{vnum,1} = vcut(leftright);
                    end
                end
                c = p;
                p = tree.Parent(p);
            end

            % Create label listing any variable with finite limits
            txt = getString(message('stats:classregtree:view:AtThisNode'));
            for j=1:size(lims,1);
                L1 = lims{j,1};
                L2 = lims{j,2};
                if ~iscell(L1) && isinf(L1) && isinf(L2)
                    continue
                end
                vnum = lims{j,3};

                if iscell(L1)
                    cats = L1{1};
                    if length(cats)==1
                        txt = sprintf('%s\n%s = %s',txt,varnames{vnum},num2str(cats,'%g '));
                    else
                        txt = sprintf('%s\n%s %s (%s)',txt,varnames{vnum},...
                            getString(message('stats:classregtree:disp:ElementInSet')),...
                            num2str(cats,'%g '));
                    end
                elseif isinf(L1)
                    txt = sprintf('%s\n%s < %g',txt,varnames{vnum},L2);
                elseif isinf(L2)
                    txt = sprintf('%s\n%g <= %s',txt,L1,varnames{vnum});
                else
                    txt = sprintf('%s\n%g <= %s < %g',txt,L1,varnames{vnum},L2);
                end
            end
        end
    else
        txt = '';
    end

    % Add label
    if ~isempty(txt)
        x = X(node);
        y = Y(node);
        xlim = get(gca,'xlim');
        ylim = get(gca,'ylim');
        if x<mean(xlim)
            halign = 'left';
            dx = 0.02;
        else
            halign = 'right';
            dx = -0.02;
        end
        if y<mean(ylim)
            valign = 'bottom';
            dy = 0.02;
        else
            valign = 'top';
            dy = -0.02;
        end
        h = text(x+dx*diff(xlim),y+dy*diff(ylim),txt,'Interpreter','none');
        yellow = [1 1 .85];
        set(h,'backgroundcolor',yellow,'margin',3,'edgecolor','k',...
            'HorizontalAlignment',halign,'VerticalAlignment',valign,...
            'tag','LinescanText','ButtonDownFcn',@startmovetips);
        line(x,y,'Color',yellow,'Marker','.','MarkerSize',20,...
            'Tag','LinescanMarker');
    end
end
end

% ----------------------------------------------
function startmovetips(varargin)
%STARTMOVETIPS Start movement of node tips.

f = gcbf;
set(f,'WindowButtonUpFcn',@donemovetips,...
      'WindowButtonMotionFcn',@showmovetips,...
      'Interruptible','off','BusyAction','queue');

o = gcbo;
p1 = get(f,'CurrentPoint');
a = get(f,'CurrentAxes');
ud = get(a,'UserData');
ud(1:2) = {o p1};
set(a,'UserData',ud);
end

% ----------------------------------------------
function showmovetips(varargin)
%SHOWMOVETIPS Show current movement of node tips.
domovetips(0,varargin{:});
end

% ----------------------------------------------
function donemovetips(varargin)
%DONEMOVETIPS Finish movement of node tips.
domovetips(1,varargin{:});
end

% ----------------------------------------------
function domovetips(alldone,varargin)
%DOMOVETIPS Carry out movement of node tips.

f = gcbf;
if alldone
   set(f,'WindowButtonUpFcn','','WindowButtonMotionFcn','');
end
a = get(f,'CurrentAxes');
ud = get(a,'UserData');
o = ud{1};
p1 = ud{2};
p2 = get(f,'CurrentPoint');
p0 = get(a,'Position');
pos = get(o,'Position');
dx = (p2(1)-p1(1)) * diff(get(a,'XLim')) / p0(3);
dy = (p2(2)-p1(2)) * diff(get(a,'YLim')) / p0(4);
pos(1) = pos(1) + dx;
pos(2) = pos(2) + dy;
set(o,'Position',pos);
ud{2} = p2;
set(a,'UserData',ud);
end

% ----------------------------------------------
function resize(varargin)
%RESIZE Resize figure showing decision tree.
layoutfig(gcbf)
end

% ----------------------------------------------
function layoutfig(f)
%LAYOUTFIG Layout figure contents

set(f,'Units','points');
fpos = get(f,'Position');

% Resize frame
h = findobj(f,'Tag','frame');
frpos = get(h,'Position');
frpos(2) = fpos(4) - frpos(4);
frpos(3) = fpos(3);
set(h,'Position',frpos);

% Resize controls inside frame
tags = {'clicktext'  'clicklist'  'magtext' 'maglist' ...
        'pruneframe' 'prunelabel' 'prunelev'};
mult = [1.6          1.35         1.6       1.35      ...
        1.7          1.6          1.6];
for j=1:length(tags)
   h = findobj(f,'Tag',tags{j});
   p = get(h,'Position');
   if j==1, theight = p(4); end
   p(2) = fpos(4) - mult(j)*theight;
   set(h,'Position',p);
end

h = findobj(f,'Tag','grow');
p = get(h,'Position');
p(2) = frpos(2)+2;
set(h,'Position',p);
h = findobj(f,'Tag','prune');
p(2) = p(2)+p(4);
set(h,'Position',p);

% Resize sliders
hh = findobj(f,'Tag','hslider');
hv = findobj(f,'Tag','vslider');
p1 = get(hh,'Position');
sw = p1(4);
p1(3) = frpos(3) - sw - 1;
set(hh,'Position',p1);
p2 = get(hv,'Position');
p2(1) = frpos(3) - sw - 1;
p2(4) = frpos(2) - sw - 1;
set(hv,'Position',p2);
if isequal(get(hh,'Visible'),'off')
   sw = 0;
end

% Resize axes
h = get(f,'CurrentAxes');
p = [0, sw, frpos(3)-sw, frpos(2)-sw];
set(h,'Position',p);
end

% ------------------------------------------
function domenu(varargin)
%DOMENU Carry out menu actions for tree viewer.

o = gcbo;
f = gcbf;
t = get(o,'Tag');
switch(t)
 % Change display from full tree to unpruned nodes or vice versa
 case {'menufull' 'menuunpr'}
   ischecked = isequal(get(o,'Checked'),'on');
   isfull = isequal(t,'menufull');
   if isfull
      dofull = ~ischecked;
   else
      dofull = ischecked;
   end
   mfull = findobj(f,'Type','uimenu','Tag','menufull');
   munpr = findobj(f,'Type','uimenu','Tag','menuunpr');
   if dofull
      set(mfull,'Checked','on');
      set(munpr,'Checked','off');
   else
      set(mfull,'Checked','off');
      set(munpr,'Checked','on');
   end
   doprunegraph(f,dofull);
   dozoom(f);

 % Turn on/off branch labels
 case 'menubr'
   curval = get(o,'Checked');
   if isequal(curval,'on')
      set(o,'Checked','off');
      h = findobj(f,'Type','text','Tag','branchlabel');
      set(h,'Visible','off');
   else
      set(o,'Checked','on');
      h = findobj(f,'Type','text','Tag','branchlabel');
      set(h,'Visible','on');
   end

 % Turn on/off leaf labels
 case 'menuleaf'
   curval = get(o,'Checked');
   if isequal(curval,'on')
      set(o,'Checked','off');
      h = findobj(f,'Type','text','Tag','leaflabel');
      set(h,'Visible','off');
   else
      set(o,'Checked','on');
      h = findobj(f,'Type','text','Tag','leaflabel');
      set(h,'Visible','on');
   end
end
end

% ------------------------------------------
function doprunegraph(f,dofull)
%DOPRUNEGRAPH Adjust graph to show full/pruned setting

a = get(f,'CurrentAxes');
h1 = findobj(a,'Type','line','Tag','prunednode');
h2 = findobj(a,'Type','line','Tag','prunedconnection');

% Figure out whether to show full tree
if nargin<2
   o = findobj(f,'Type','uimenu','Tag','menufull');
   dofull = isequal(get(o,'Checked'),'on');
end

% Adjust graph
if dofull
   set(h1,'Visible','on');
   set(h2,'Visible','on');
   xlim = get(a,'XLim');
   ylim = get(a,'YLim');
   bigxlim = 0:1;
   bigylim = 0:1;
else
   set(h1,'Visible','off');
   set(h2,'Visible','off');
   h1 = findobj(f,'Type','line','Tag','leaf');
   h2 = findobj(f,'Type','line','Tag','branch');
   x1 = get(h1,'XData');
   y1 = get(h1,'YData');
   y2 = get(h2,'YData');
   dx = 1 / (1+length(x1));
   ally = sort(unique([y1(:); y2(:)]));
   if length(ally)>1
      dy = 0.5 * (ally(2)-ally(1));
   else
      dy = 1-ally;
   end
   xlim = [min(x1)-dx, max(x1)+dx];
   ylim = [min(ally)-dy, max(ally)+dy];
   bigxlim = 0:1;
   bigylim = [ylim(1) 1];
end
axislistener(a,false);
set(a,'XLim',xlim,'YLim',ylim);
axislistener(a,true);
hh = findobj(f,'Tag','hslider');
set(hh,'UserData',bigxlim);
hv = findobj(f,'Tag','vslider');
set(hv,'UserData',bigylim);
end

% ------------------------------------------
function domagnif(varargin)
%DOMAGNIF React to magnification level change

f = gcbf;
o = gcbo;

% We need sliders if the magnification level is > 100%
h = [findobj(f,'Tag','hslider'), findobj(f,'Tag','vslider')];
maglevel = get(o,'Value');
if maglevel==1
   set(h,'Visible','off');
else
   set(h,'Visible','on');
end

% Adjust layout if necessary
resize;

% Adjust axes contents
dozoom(f);

% Remove custom zoom amount from list if not in use
if maglevel<=4
   adjustcustomzoom(o,false);
end

% Turn off any manual zooming
zoom(f,'off');
end

% --------------------------------------------------
function adjustcustomzoom(o,add)
%ADJUSTCUSTOMZOOM Add or remove special custom magnification level
nchoices = size(get(o,'String'),1);
choices = '100%|200%|400%|800%';
if ~add && nchoices~=4
   set(o,'String',choices);
elseif add && nchoices~=5
   choices = [choices '|' 'Custom'];
   set(o,'String',choices);
end
end

% ------------------------------------------
function dozoom(f)
%DOZOOM Adjust axes contents to match magnification settings

a = get(f,'CurrentAxes');
hh = findobj(f,'Tag','hslider');
hv = findobj(f,'Tag','vslider');
hm = findobj(f,'Tag','maglist');

% Get information about x/y ranges and current midpoint
bigxlim = get(hh,'UserData');
bigylim = get(hv,'UserData');
xlim = get(a,'XLim');
ylim = get(a,'YLim');
currx = (xlim(1)+xlim(2))/2;
curry = (ylim(1)+ylim(2))/2;

% How much are we magnifying each axis?
magfact = [1 2 4 8];
mag = get(hm,'Value');
if mag<=4
   magfact = magfact(mag)*ones(1,2);
else
   magfact = [diff(bigxlim)/diff(xlim), diff(bigylim)/diff(ylim)];
end
magfact = max(magfact,1);

if all(magfact==1)                 % no magnification
   xlim = bigxlim;
   ylim = bigylim;
else                               % magnify by showing a subset of range
   magfact = max(magfact,1.01);
   dx = diff(bigxlim)/magfact(1);
   dy = diff(bigylim)/magfact(2);
   xval = max(bigxlim(1), min(bigxlim(2)-dx, currx-dx/2));
   xlim = xval + [0,dx];
   yval = max(bigylim(1), min(bigylim(2)-dy, curry-dy/2));
   ylim = yval + [0,dy];
   set(hh,'Min',bigxlim(1),'Max',bigxlim(2)-dx,'Value',xval);
   set(hv,'Min',bigylim(1),'Max',bigylim(2)-dy,'Value',yval);
end
axislistener(a,false);
set(a,'XLim',xlim,'YLim',ylim);
axislistener(a,true);
end

% ------------------------------------------
function dopan(varargin)
%DOPAN Pan around magnified tree display

f = gcbf;
a = get(f,'CurrentAxes');
o = gcbo;
val = get(o,'Value');

axislistener(a,false);
if isequal(get(o,'Tag'),'hslider')
   xlim = get(a,'XLim');
   xlim = xlim + (val-xlim(1));
   set(a,'XLim',xlim);
else
   ylim = get(a,'YLim');
   ylim = ylim + (val-ylim(1));
   set(a,'YLim',ylim);
end
axislistener(a,true);
end

% -----------------------------------------
function axislistener(a,enable)
%AXISLISTENER Enable or disable listening to axis limit changes

f = get(a,'Parent');
ud = get(a,'UserData');
if enable
   % Start listening to axis limit changes
   list1 = addlistener(a, 'XLim', 'PostSet', @(src,evt) customzoom(f));
   list2 = addlistener(a, 'YLim', 'PostSet', @(src,evt) customzoom(f));
   ud(3:4) = {list1 list2};
else
   % Delete listeners
   for j=3:4
      lstnr = ud{j};
      if ~isempty(lstnr), delete(lstnr); end
   end
   ud(3:4) = {[]};
end
set(a,'UserData',ud);
end

% -----------------------------------------
function customzoom(f)
%CUSTOMPAN Deal with panning of a zoomed tree view

a = get(f,'CurrentAxes');
xlim = get(a,'XLim');
ylim = get(a,'YLim');

hh = findobj(f,'Tag','hslider');
hv = findobj(f,'Tag','vslider');
hm = findobj(f,'Tag','maglist');

bigxlim = get(hh,'UserData');
bigylim = get(hv,'UserData');
magfact = [1 2 4 8];

% Figure out if we have a standard zoom amount (100%, 200%, etc.) or
% a custom zoom amount
xratio = diff(bigxlim) / diff(xlim);
yratio = diff(bigylim) / diff(ylim);
standard = abs(xratio-yratio)<=0.02 & abs(xratio-round(xratio))<=0.02 ...
                                    & abs(yratio-round(yratio))<=0.02;
if standard
   xratio = round(xratio);
   standard = ismember(xratio,magfact);
end

% Update the magnification setting
if standard
   set(hm,'Value',find(magfact==xratio));
   adjustcustomzoom(hm,false);
   if xratio==1
      h = [findobj(f,'Tag','hslider'), findobj(f,'Tag','vslider')];
      set(h,'Visible','off');
   end
else
   adjustcustomzoom(hm,true);
   set(hm,'Value',5);
   h = [findobj(f,'Tag','hslider'), findobj(f,'Tag','vslider')];
   set(h,'Visible','on');
end

dozoom(f);
end
