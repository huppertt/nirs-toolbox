classdef TreeImpl
    
%   Copyright 2012-2014 The MathWorks, Inc.

    properties
        D = []; % number of predictors
        Children = [];
        ClassCount = [];
        ClassNames = {}; % empty unless made from classregtree
        ClassProb = [];
        CutCategories
        CutPoint = [];
        CutVar = [];
        HasUnsplit = [];
        IsBranch = [];
        NodeMean = [];
        NodeProb = [];
        NodeRisk = [];
        NodeSize = [];
        Parent = [];
        PruneList = [];
        PruneAlpha = [];
        SplitGain = [];
        SurrCutCategories = [];
        SurrCutFlip = [];
        SurrCutPoint = [];
        SurrCutVar = [];
        SurrSplitGain = [];
        SurrVarAssoc = [];
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        CatSplit;
        CutType;
        SurrCutType;
    end
    
    methods
        function a = get.CatSplit(this)
            branches = this.IsBranch;
            notNumericCut = isnan(this.CutPoint);
            catCut = branches & notNumericCut;
            if any(catCut)
                a = this.CutCategories(catCut,:);
            else
                a = {};
            end
        end
        
        function a = get.CutType(this)
            branches = this.IsBranch;
            numericCut = ~isnan(this.CutPoint);
            N = length(branches);
            a = repmat({''},N,1);
            a(branches & numericCut)  = {'continuous'};
            a(branches & ~numericCut) = {'categorical'};
        end
        
        function a = get.SurrCutType(this)
            if isempty(this.SurrCutPoint)
                a = {};
            else
                a = repmat({{}},size(this.IsBranch,1),1);
                branches = find(this.IsBranch);
                for b=1:numel(branches)
                    node = branches(b);
                    cutpoint = this.SurrCutPoint{node};
                    if ~isempty(cutpoint)
                        nodecut = repmat({''},1,numel(cutpoint));
                        numericCut = ~isnan(cutpoint);
                        nodecut(numericCut) = {'continuous'};
                        nodecut(~numericCut) = {'categorical'};
                        a{node} = nodecut;
                    end
                end
            end
        end
    end
    
    methods(Access=protected)
        function this = TreeImpl()
        end

        function this = pruneNodes(this,branches)
            % Prune away children of these nodes
            N = size(this.Children,1);
            
            % Find children of these branches and remove them
            parents = branches;
            tokeep = true(N,1);
            kids = [];
            while(true)
                newkids = this.Children(parents,:);
                newkids = newkids(:);
                newkids = newkids(newkids>0 & ~ismember(newkids,kids));
                if isempty(newkids)
                    break;
                end
                kids = [kids; newkids];
                tokeep(newkids) = false;
                parents = newkids;
            end
            
            % Convert branches to leaves by removing split rule and children
            this.CutVar(branches) = 0;
            this.CutPoint(branches) = NaN;
            this.CutCategories(branches,:) = {[]};
            this.Children(branches,:) = 0;
            this.HasUnsplit(branches) = false;
            this.IsBranch(branches) = false;
            this.SplitGain(branches) = 0;

            if ~isempty(this.SurrCutVar)
                this.SurrCutVar(branches) = {[]};
                this.SurrCutPoint(branches) = {[]};
                this.SurrCutCategories(branches) = {{}};
                this.SurrCutFlip(branches) = {[]};
                this.SurrSplitGain(branches) = {[]};
                this.SurrVarAssoc(branches) = {[]};
            end
            
            % Get new node numbers from old node numbers
            ntokeep = sum(tokeep);
            nodenums = zeros(N,1);
            nodenums(tokeep) = (1:ntokeep)';
            
            % Nodes to remove
            remove = ~tokeep;
            
            % Update node numbers
            this.Parent(remove)                = [];
            this.Children(remove,:)            = [];
            mask = this.Parent>0;
            this.Parent(mask) = nodenums(this.Parent(mask));
            mask = this.Children>0;
            this.Children(mask) = nodenums(this.Children(mask));
            
            % Get rid of nodes
            this.ClassCount(remove,:)          = [];
            this.ClassProb(remove,:)           = [];
            this.CutCategories(remove,:)       = [];
            this.CutPoint(remove)              = [];
            this.CutVar(remove)                = [];
            this.HasUnsplit(remove)            = [];
            this.IsBranch(remove)              = [];
            this.NodeMean(remove)              = [];
            this.NodeProb(remove)              = [];
            this.NodeRisk(remove)              = [];
            this.NodeSize(remove)              = [];
            this.SplitGain(remove)             = [];
            
            if ~isempty(this.SurrCutVar)
                this.SurrCutVar(remove)        = [];
                this.SurrCutPoint(remove)      = [];
                this.SurrCutCategories(remove) = [];
                this.SurrCutFlip(remove)       = [];
                this.SurrSplitGain(remove)     = [];
                this.SurrVarAssoc(remove)      = [];
            end
        end
    end
        
    methods
        function this = prune(this,varargin)
            verbose = 0;
            
            args = {'forceprune' 'criterion' 'cost' 'level' 'nodes' 'alpha'};
            defs = {       false          ''     []      []      []      []};
            [force,crit,cost,level,nodes,alpha] = ...
                internal.stats.parseArgs(args,defs,varargin{:});

            dolevel = ~isempty(level);
            donodes = ~isempty(nodes);
            doalpha = ~isempty(alpha);
            
            if sum(dolevel+donodes+doalpha)>1
                error(message('stats:classreg:learning:impl:TreeImpl:prune:TooManyPruningOptions'));
            end
            
            doprune = any([dolevel donodes doalpha]);
                
            if ~isempty(crit) && (~ischar(crit) || ~isvector(crit))
                error(message('stats:classreg:learning:impl:TreeImpl:prune:BadCrit'));
            end
                        
            if dolevel
                if ~isnumeric(level) || ~isscalar(level) || level<0
                    error(message('stats:classreg:learning:impl:TreeImpl:prune:BadLevel'));
                end
                level = ceil(level);
            end

            if doalpha
                if ~isnumeric(alpha) || ~isscalar(alpha) || any(alpha<0)
                    error(message('stats:classreg:learning:impl:TreeImpl:prune:BadAlpha'));
                end
            end
               
            if donodes
                if ~isnumeric(nodes) || ~isvector(nodes) || any(nodes<1)
                    error(message('stats:classreg:learning:impl:TreeImpl:prune:BadNodes'));
                end
                nodes = ceil(nodes);
            end
            
            if isempty(this.PruneList) || isempty(this.PruneAlpha) || force
                [prunelist,prunealpha] = ...
                    classreg.learning.treeutils.computePruneInfo(...
                    this.ClassProb',cost,...
                    this.NodeProb,this.NodeRisk,this.HasUnsplit,...
                    this.Children',crit,verbose);
                this.PruneList         = prunelist;
                this.PruneAlpha        = prunealpha;
            end

            % If no pruning, simply compute the pruning sequence
            if ~doprune
                return;
            end
            
            if dolevel && level>max(this.PruneList)
                warning(message('stats:classreg:learning:impl:TreeImpl:prune:LevelTooLarge', ...
                    level, max(this.PruneList)));
            end
            
            if doalpha && alpha>max(this.PruneAlpha)
                warning(message('stats:classreg:learning:impl:TreeImpl:prune:AlphaTooLarge', ...
                    sprintf('%g',alpha), sprintf('%g',max(this.PruneAlpha))));
            end
            
            if doalpha
                level = find(this.PruneAlpha<=alpha,1,'last') - 1;
            end
            
            % Find nodes whose children need to be pruned
            if ~isempty(level)
                nodes = find(this.IsBranch & this.PruneList<=level);
            end
            
            % Prune nodes
            pruned = false;
            if ~isempty(nodes)
                this = pruneNodes(this,nodes);
                pruned = true;
            end

            % Recompute pruning sequence
            if pruned
                [prunelist,prunealpha] = ...
                    classreg.learning.treeutils.computePruneInfo(...
                    this.ClassProb',cost,...
                    this.NodeProb,this.NodeRisk,this.HasUnsplit,...
                    this.Children',crit,verbose);
                this.PruneList         = prunelist;
                this.PruneAlpha        = prunealpha;
            end
        end
        
        function subtrees = processSubtrees(this,subtrees)
            if ~strcmpi(subtrees,'all') && ...
                    (~isnumeric(subtrees) || ~isvector(subtrees) ...
                    || any(subtrees<0) || any(diff(subtrees)<0))
                error(message('stats:classreg:learning:impl:TreeImpl:processSubtrees:BadSubtrees'));
            end
            if isscalar(subtrees) && subtrees==0
                return;
            end
            prunelevs = this.PruneList;
            if isempty(prunelevs)
                error(message('stats:classreg:learning:impl:TreeImpl:processSubtrees:NoPruningInfo'));
            end
            if ischar(subtrees)
                subtrees = min(prunelevs):max(prunelevs);
            end
            subtrees = ceil(subtrees);
            if subtrees(end)>max(prunelevs)
                error(message('stats:classreg:learning:impl:TreeImpl:processSubtrees:SubtreesTooBig'));
            end
        end
                
        function tree = findsubtree(this,alpha0)
            adjfactor = 1 + 100*eps;
            alpha = this.PruneAlpha;
            tree = zeros(size(alpha0));
            for j=1:length(alpha0);
                tree(j) = sum(alpha <= alpha0(j)*adjfactor);
            end
            tree = tree - 1;
        end
        
        function nleaf = countLeaves(this,subtrees)
            N = numel(this.PruneAlpha);
            if N==0
                nleaf = sum(~this.IsBranch);
            else
                if strcmp(subtrees,'all')
                    subtrees = 0:N-1;
                end
                N = numel(subtrees);
                nleaf = zeros(N,1);
                for n=1:N
                    t = prune(this,'level',subtrees(n));
                    nleaf(n) = sum(~t.IsBranch);
                end
            end
        end
            
        function n = findNode(this,X,catpred,subtrees)
            verbose = 0;
            
            if ~isfloat(X) || ~ismatrix(X)
                error(message('stats:classreg:learning:impl:TreeImpl:findNode:BadX'));
            end
            internal.stats.checkSupportedNumeric('X',X);
            
            p = size(X,2);            
            if p~=this.D
                error(message('stats:classreg:learning:impl:TreeImpl:findNode:BadXSize', this.D));
            end
        
            iscat = false(size(X,2),1);
            iscat(catpred) = true;
            
            % Validate subtrees
            subtrees = processSubtrees(this,subtrees);
                        
            n = classreg.learning.treeutils.findNode(X,...
                subtrees,this.PruneList,...
                this.Children',iscat,...
                this.CutVar,this.CutPoint,this.CutCategories,...
                this.SurrCutFlip,this.SurrCutVar,...
                this.SurrCutPoint,this.SurrCutCategories,...
                verbose);
        end
        
        
        function [imp,nsplit] = predictorImportance(this,varargin)
            imp = zeros(1,this.D);
            nsplit = zeros(1,this.D);

            nBranch = sum(this.IsBranch);
            if nBranch==0
                return;
            end
            
            splitGain = this.SplitGain;
            cutVar = this.CutVar;
            for d=1:this.D
                imp(d) = sum(splitGain(cutVar==d));
                nsplit(d) = sum(cutVar==d);
            end
            
            surrCutVar = this.SurrCutVar;
            surrSplitGain = this.SurrSplitGain;
            if ~isempty(surrCutVar)
                M = numel(surrCutVar);
                for m=1:M
                    thisCutVar = surrCutVar{m};
                    thisSplitGain = surrSplitGain{m};
                    imp(thisCutVar) = imp(thisCutVar) + thisSplitGain;
                    nsplit(thisCutVar) = nsplit(thisCutVar) + 1;
                end
            end
            
            imp = imp/nBranch;
        end    
                
        function ma = meanSurrVarAssoc(this,j)
            if nargin>1 
                validateNodes(this,j);
            end
            
            if nargin<2
                j = 1:numel(this.SurrCutVar);
            end
            
            % Keep only branch nodes
            isbr = this.IsBranch(j);
            j = j(isbr);
            
            % Init
            N = numel(j);
            p = this.D;
            ma = zeros(p);
            nsplit = zeros(p,1);
            
            % Get the association matrix and lists of best and surrogate predictors
            a = this.SurrVarAssoc(j);
            bestvar = this.CutVar(j);
            surrvar = this.SurrCutVar(j);
            
            % Loop over optimal splits. Increase the split count by 1 for every node.
            for i=1:N
                n = bestvar(i);
                nsplit(n) = nsplit(n) + 1;
                m = surrvar{i};
                if ~isempty(m)
                    ma(n,m) = ma(n,m) + a{i};
                end
            end
            
            % Divide cumulative association by the number of optimal splits on each
            % predictor
            gt0 = nsplit>0;
            if any(gt0)
                ma(gt0,:) = bsxfun(@rdivide,ma(gt0,:),nsplit(gt0));
            end
            
            % Make sure the diagonal elements are 1
            ma(1:p+1:end) = 1;
        end
    
        function validateNodes(this,j)
            numnodes = size(this.Children,1);
            if islogical(j)
                ok = (numel(j) <= numnodes);
            else
                ok = ismember(j,1:numnodes);
                ok = all(ok(:));
            end
            if ~ok
                error(message('stats:classreg:learning:impl:TreeImpl:validateNodes:BadNodes', numnodes, numnodes));
            end
        end
            
        function view(this,classnames,nodevalue,varnames,htmlhelp,varargin)
            args = {'mode' 'prunelevel'};
            defs = {'text'           []};
            [mode,prunelevel] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            if     strncmpi(mode,'text',length(mode))
                viewText(this,classnames,nodevalue,varnames,prunelevel);
            elseif strncmpi(mode,'graph',length(mode))
                viewGraph(this,classnames,nodevalue,varnames,prunelevel,htmlhelp);
            else
                error(message('stats:classreg:learning:impl:TreeImpl:view:BadViewMode'));
            end
        end        
        
        function viewText(this,classnames,nodevalue,varnames,prunelevel)
            
            if isempty(prunelevel)
                t = this;
            else
                t = prune(this,'level',prunelevel);
            end
            
            if ~isempty(classnames)
                nodevalue = classnames(nodevalue);
            end
            
            isLoose = strcmp(get(0,'FormatSpacing'),'loose');
            if (isLoose), fprintf('\n'); end
            
            % Get some information about the whole tree
            maxnode = size(t.Children,1);
            nd = 1 + floor(log10(maxnode)); % number of digits for node number
            isregression = isempty(classnames);
            if isregression
                fprintf(getString(message('stats:classregtree:disp:DecisionTreeForRegression')));
            else
                fprintf(getString(message('stats:classregtree:disp:DecisionTreeForClassification')));
            end
            
            % Display information about each node. nodevalue is cellstr for
            % classification tree and numeric vector for regression tree.
            for j=1:maxnode
                if any(t.Children(j,:))
                    % branch node
                    vnum = t.CutVar(j);
                    vtype = t.CutType{j};
                    vname = varnames{vnum};
                    kids = t.Children(j,:);
                    if isempty(classnames)
                        Yfit = nodevalue(j);
                        Yfit = num2str(Yfit,'%g');
                    else
                        Yfit = nodevalue{j};
                    end
                    if strcmp(vtype,'continuous')
                        cut = t.CutPoint(j);
                        condleft = sprintf('%s<%g',vname,cut);
                        condright = sprintf('%s>=%g',vname,cut);
                        fprintf('%*d  %s\n',nd,j,getString(message('stats:classregtree:disp:TreeBranch',...
                            condleft,kids(1),condright,kids(2),Yfit)));
                    else
                        cut = t.CutCategories(j,:);
                        cats = cut{1};
                        if isscalar(cats)
                            condleft = sprintf('%s=%g',vname,cats);
                        else
                            set = deblank(num2str(cats,'%g '));
                            condleft = sprintf('%s %s {%s}',vname,getString(message('stats:classregtree:disp:ElementInSet')),set);
                        end
                        cats = cut{2};
                        if isscalar(cats)
                            condright = sprintf('%s=%g',vname,cats);
                        else
                            set = deblank(num2str(cats,'%g '));
                            condright = sprintf('%s %s {%s}',vname,getString(message('stats:classregtree:disp:ElementInSet')),set);
                        end
                        fprintf('%*d  %s\n',nd,j,getString(message('stats:classregtree:disp:TreeBranch',...
                            condleft,kids(1),condright,kids(2),Yfit)));
                    end
                else
                    % terminal node, display fit (regression) or class assignment
                    if isregression
                        fprintf(sprintf('%s  %s %s\n','%*d',getString(message('stats:classregtree:disp:FittedResponse')),'%g'),nd,j,nodevalue(j));
                    else
                        fprintf(sprintf('%s  %s %s\n','%*d',getString(message('stats:classregtree:disp:PredictedClass')),'%s'),nd,j,nodevalue{j});
                    end
                end
            end
            if (isLoose), fprintf('\n'); end
            
        end
        
        function outfig = viewGraph(this,classnames,nodevalue,varnames,curlevel,htmlhelp)
            
            if isempty(curlevel)
                curlevel = 0;
            end
            
            doclass = ~isempty(classnames);
            
            function helpviewer(varargin)
                helpview([docroot htmlhelp]);
            end
            
            % Create empty figure and axes to receive tree display
            fig = classreg.learning.treeutils.TreeDrawer.setupfigure(doclass);
            try
                classreg.learning.treeutils.TreeDrawer.adjustmenu(fig,@helpviewer);
            catch me
                error(message('stats:classreg:learning:impl:TreeImpl:view:AdjustMenuFails', me.message));
            end
            
            % Save the full tree
            fulltree = this;
            
            % Draw tree
            [X,Y] = classreg.learning.treeutils.TreeDrawer.drawtree(...
                this,fig,nodevalue,varnames,curlevel,classnames);
            
            % Save information for call-backs
            set(fig,'ButtonDownFcn',@classreg.learning.treeutils.TreeDrawer.removelabels, ...
                'UserData',{X Y 0 varnames nodevalue fulltree curlevel classnames});
            
            % Update ui elements
            classreg.learning.treeutils.TreeDrawer.updateenable(fig);
            classreg.learning.treeutils.TreeDrawer.updatelevel(fig,curlevel,fulltree);
            
            if nargout>0
                outfig = fig;
            end
        end
    end
    
    methods(Static)
        function this = makeFromData(X,Y,W,useObs,doclass,catpred,splitcrit,...
                minleaf,minparent,maxsplits,nvartosample,nsurrsplit,...
                maxcat,algcat,cost,reltol,rsh)
            
            % Setting verbose to a positive scalar switches to a
            % single-thread execution for computing optimal decision
            % splits. Multithreading is preserved for sorting and NaN
            % removal.
            verbose = 0;
            
            iscat = false(size(X,2),1);
            iscat(catpred) = true;
            
            if strcmpi(nsurrsplit,'all')
                nsurrsplit = size(X,2);
            end
            
            if strcmpi(nvartosample,'all')
                nvartosample = size(X,2);
            end
            
            [~,algcat] = ismember(algcat,{'auto' 'Exact' 'PullLeft' 'PCA' 'OVAbyClass'});
            algcat = algcat-1;
            
            if doclass
                Y = Y-1; % group numbers starting at zero
            end
            
            [children,classcount,classprob,...
                cutcategories,cutpoint,cutvar,hasunsplit,isbranch,...
                nodemean,nodeprob,nodesize,parent,noderisk,splitgain,...
                surrcutcategories,surrcutflip,...
                surrcutpoint,surrcutvar,surrsplitgain,surrvarassoc] = ...
                classreg.learning.treeutils.growTree(...
                X,Y,W,useObs-1,iscat,splitcrit,...
                minleaf,minparent,maxsplits,...
                nvartosample,nsurrsplit,maxcat,algcat,cost,reltol,...
                rsh,verbose);
            
            this                   = classreg.learning.impl.TreeImpl();
            this.D                 = size(X,2);  
            this.Children          = children';
            this.ClassCount        = classcount';
            this.ClassProb         = classprob';
            this.CutCategories     = cutcategories;
            this.CutPoint          = cutpoint;
            this.CutVar            = cutvar;
            this.HasUnsplit        = hasunsplit;
            this.IsBranch          = isbranch;
            this.NodeMean          = nodemean;
            this.NodeProb          = nodeprob;
            this.NodeRisk          = noderisk;            
            this.NodeSize          = nodesize;
            this.Parent            = parent;
            this.SplitGain         = splitgain;
            this.SurrCutCategories = surrcutcategories;
            this.SurrCutFlip       = surrcutflip;
            this.SurrCutPoint      = surrcutpoint;
            this.SurrCutVar        = surrcutvar;
            this.SurrSplitGain     = surrsplitgain;
            this.SurrVarAssoc      = surrvarassoc;
        end
        
        
        function this = makeFromClassregtree(tree,calledFromLoadobj)
            if nargin<2
                calledFromLoadobj = false;
            end
            
            if ~isempty(tree.impurity)
                crit = 'impurity';
            else
                crit = 'error';
            end
            
            N = numnodes(tree);
            
            this                   = classreg.learning.impl.TreeImpl();            
            this.D                 = tree.npred;
            this.Children          = children(tree);
            
            if strcmp(type(tree),'classification')
                this.ClassCount    = classcount(tree);
            else
                this.ClassCount    = nodesize(tree);
            end
            
            this.ClassNames        = classreg.learning.internal.ClassLabel(classname(tree));
            
            if strcmp(type(tree),'classification')
                this.ClassProb     = classprob(tree);
            else
                this.ClassProb     = ones(N,1);
            end
            
            this.CutCategories     = cutcategories(tree);
            this.CutPoint          = cutpoint(tree);
            this.CutVar            = abs(tree.var);
            this.HasUnsplit        = true(N,1);
            this.IsBranch          = isbranch(tree);
            
            if strcmp(type(tree),'regression')
                this.NodeMean      = nodemean(tree);
            else
                this.NodeMean      = NaN(N,1);
            end
            
            this.NodeProb          = nodeprob(tree);
            this.NodeRisk          = risk(tree,1:N,'criterion',crit);
            this.NodeSize          = nodesize(tree);
            this.Parent            = parent(tree);
            
            this.PruneList         = tree.prunelist;
            this.PruneAlpha        = tree.alpha;
            
            r     = risk(tree,1:N,'criterion',crit);
            rdiff = risk(tree,1:N,'criterion',crit,'mode','diff');
            this.SplitGain         = rdiff;
            hasKids = find(all(children(tree)>0,2));
            
            this.SplitGain(hasKids)= rdiff(hasKids) ...
                - sum(reshape(r(children(tree,hasKids)),numel(hasKids),2),2);
            
            this.SurrCutCategories = surrcutcategories(tree);
            this.SurrCutFlip       = surrcutflip(tree);
            this.SurrCutPoint      = surrcutpoint(tree);
            this.SurrCutVar        = cellfun(@abs,tree.surrvar,'UniformOutput',false);
            
            M = numel(this.SurrCutVar);
            this.SurrSplitGain = cell(M,1);
            if M>0 && calledFromLoadobj
                warning(message('stats:classreg:learning:impl:TreeImpl:makeFromClassregtree:Pre13aTreeWithSurrogateSplits'));
            end
            for m=1:M
                this.SurrSplitGain{m} = zeros(1,numel(this.SurrCutVar{m}));
            end

            this.SurrVarAssoc      = surrvarassoc(tree);
        end
    end
    
end
