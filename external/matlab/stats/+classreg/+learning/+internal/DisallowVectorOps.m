classdef DisallowVectorOps

%   Copyright 2010 The MathWorks, Inc.


    methods(Access=protected)
        % Only derived classes can construct objects of this class.
        function this = DisallowVectorOps()
        end
    end
    
    methods(Hidden)
        function [varargout] = subsref(this,s)
            className = class(this);
            if     strcmp(s(1).type,'()')
                error(message('stats:classreg:learning:internal:DisallowVectorOps:subsref:SubscriptReferenceNotAllowed',className));
            elseif strcmp(s(1).type,'{}')
                error(message('stats:classreg:learning:internal:DisallowVectorOps:subsref:CellReferenceNotAllowed',className));
            else
                % Return default subsref to this object
                [varargout{1:nargout}] = builtin('subsref',this,s);
            end
        end
        
        function [varargout] = subsasgn(this,s,data)
            className = class(this);
            if     strcmp(s(1).type,'()')
                error(message('stats:classreg:learning:internal:DisallowVectorOps:subsasgn:SubscriptAssignmentNotAllowed',className));
            elseif strcmp(s(1).type,'{}')
                error(message('stats:classreg:learning:internal:DisallowVectorOps:subsasgn:CellAssignmentNotAllowed',className));
            else
                % Return default subsasgn to this object
                [varargout{1:nargout}] = builtin('subsasgn',this,s,data);
            end
        end
        
        function a = lt(this,varargin),         throwUndefinedError(this,'lt'); end
        function a = le(this,varargin),         throwUndefinedError(this,'le'); end
        function a = ge(this,varargin),         throwUndefinedError(this,'ge'); end
        function a = gt(this,varargin),         throwUndefinedError(this,'gt'); end
        function a = ctranspose(this,varargin), throwUndefinedError(this,'ctranspose'); end
        function a = transpose(this,varargin),  throwUndefinedError(this,'transpose'); end
        function a = permute(this,varargin),    throwUndefinedError(this,'permute'); end
        function a = reshape(this,varargin),    throwUndefinedError(this,'reshape'); end
        function [a,b] = sort(this,varargin),   throwUndefinedError(this,'sort'); end
        
        function a = cat(this,varargin),        throwNoCatError(this); end
        function a = horzcat(this,varargin),    throwNoCatError(this); end
        function a = vertcat(this,varargin),    throwNoCatError(this); end
        function a = repmat(this,varargin),     throwNoCatError(this); end

        function throwUndefinedError(this,s)
            error(message('stats:classreg:learning:internal:DisallowVectorOps:UndefinedMethod', s, class( this )));
        end
        
        function throwNoCatError(this)
            error(message('stats:classreg:learning:internal:DisallowVectorOps:ConcatenationNotAllowed',class(this)));
        end
    end
    
end

