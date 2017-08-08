classdef CompactTreeImpl
%CompactTreeImpl Decision tree implementation. This class is used for
% pre-13a implementations based on classregtree. If you delete it, you need
% to change loadobj implementations for the compact tree objects.

%   Copyright 2010 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected)
        Tree = [];
    end
    
    methods
        function this = CompactTreeImpl(tree)
            this.Tree = tree;
        end

        function view(this,varargin)
            mode = internal.stats.parseArgs({'mode'},{'text'},varargin{:});
            if     strncmpi(mode,'text',length(mode))
                disp(this.Tree);
            elseif strncmpi(mode,'graph',length(mode))
                view(this.Tree);
            else
                error(message('stats:classreg:learning:impl:CompactTreeImpl:view:BadViewMode'));
            end
        end
        
        function [cost,secost,nleaf,bestlevel] = ...
                loss(this,X,Y,mode,subtrees,treesize,varargin)
            if size(X,1)~=numel(classreg.learning.internal.ClassLabel(Y))
                error(message('stats:classreg:learning:impl:CompactTreeImpl:loss:SizeXYMismatch'));
            end
            
            % Check treesize argument
            if ~ischar(treesize) || ~(treesize(1)=='s' || treesize(1)=='m')
                error(message('stats:classreg:learning:impl:CompactTreeImpl:loss:BadTreeSize'));
            end
            
            % Get info for pruned trees, then determine bestlevel.
            [cost,secost,nleaf] = test(this.Tree,mode,X,Y,varargin{:});
            if ~ischar(subtrees)
                cost = cost(1+subtrees);
                secost = secost(1+subtrees);
                nleaf = nleaf(1+subtrees);
            end

            % Find bestlevel
            if nargout>3
                [mincost,minloc] = min(cost);
                if isequal(treesize(1),'m')
                    cutoff = mincost * (1 + 100*eps);
                else
                    cutoff = mincost + secost(minloc);
                end
                bestlevel = subtrees(find(cost<=cutoff,1,'last'));
            end
        end
        
        function subtrees = processSubtrees(this,subtrees)
            % Check subtrees
            if ~strcmpi(subtrees,'all') ...
                    && (~isnumeric(subtrees) || ~isvector(subtrees) || any(subtrees<0))
                error(message('stats:classreg:learning:impl:CompactTreeImpl:processSubtrees:BadSubtrees'));
            end
            if isscalar(subtrees) && subtrees==0
                return;
            end
            prunelevs = prunelist(this.Tree);            
            if isempty(prunelevs)
                error(message('stats:classreg:learning:impl:CompactTreeImpl:processSubtrees:NoPruningInfo'));
            end
            if ischar(subtrees)
                subtrees = min(prunelevs):max(prunelevs);
            end
            subtrees = ceil(subtrees);
            if any(subtrees>max(prunelevs))
                error(message('stats:classreg:learning:impl:CompactTreeImpl:processSubtrees:SubtreesTooBig'));
            end
        end
    end

end
