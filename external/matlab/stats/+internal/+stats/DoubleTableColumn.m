classdef DoubleTableColumn < double
% DoubleTableColumn
%     The DoubleTableColumn class represents a column vector of doubles
%     to be used in a table display in a dataset. It provides a method
%     to indicate that some values should be considered "absent" and not
%     displayed in the table. An example might be an anova table, where
%     there is a column of F statistics but the "Error" row should not
%     have an F statistic.
%
%     See also DoubleTableColumn.DoubleTableColumn.
    

%   Copyright 2011-2013 The MathWorks, Inc.

    properties
        absent
    end
    
    methods(Access='public')
        function d = DoubleTableColumn(data,absent)
% DoubleTableColumn.DoubleTableColumn Constructor
%     C = DoubleTableColumn(VALS,ABSENT) accepts a column vector VAL
%     of double values and a logical vector ABSENT of the same size
%     as V. The result C displays blanks in places where ABSENT has
%     the value TRUE.
%
%     The output C is suitable for use as a variable in a dataset
%     when some variable values should not be present.

            if nargin==0
                data = [];
            end
            
            if ~isfloat(data) || ~iscolumn(data)
                error(message('stats:internal:DoubleTableColumn:BadData'));
            end
            if nargin<2 || isempty(absent)
                absent = false(size(data));
            elseif ~isequal(size(data),size(absent)) || ~islogical(absent)
                error(message('stats:internal:DoubleTableColumn:BadAbsent'));
            end
            d = d@double(data);
            d.absent = absent;
        end
        
        function disp(d)
            varStr = evalc('disp(double(d))');
            
            % varStr is a single row with \n delimiting the chars for
            % each row of var.
            loc = [0 find(varStr==10)];
            [n,~] = size(d); % already checked is 2D
            % Split the \n-delimited string into a char matrix.
            len = diff(loc);
            varChars = repmat(' ',size(d,1),max(len)-1);
            for i = 1:n
                celChars = varStr(loc(i)+1:loc(i+1)-1);
                if ~isempty(celChars) % avoid 0x0 coming from strtrim
                    varChars(i,1:length(celChars)) = celChars;
                end
            end
            varChars(d.absent,:) = ' ';
            disp(varChars)
        end
        
        function a = num2str(d,varargin)
            a = num2str(double(d),varargin{:});
            a(d.absent,:) = ' ';
        end
        
        function sref = subsref(obj,s)
            switch s(1).type
                case '.'
                    switch s(1).subs
                        case 'data'
                            sref = double(obj);
                        case 'absent'
                            sref = obj.absent;
                            if length(s)>1 && strcmp(s(2).type, '()')
                                sref = subsref(sref,s(2:end));
                            end
                        otherwise
                            error(message('stats:internal:DoubleTableColumn:BadField'));
                    end
                case '()'
                    d = double(obj);
                    if ~isscalar(s)
                        error(message('stats:internal:DoubleTableColumn:BadParenIndex'));
                    end
                    sref = subsref(d,s);
% *** it may be better to return a double
%                     sref = DoubleTableColumn(subsref(d,s),...
%                                     subsref(obj.absent,s));
                otherwise
                    error(message('stats:internal:DoubleTableColumn:BadCellIndex'));
            end
        end
        function idx = subsidx(s)
            idx = double(s);
        end
        
        function newobj = horzcat(varargin)
            c1 = cellfun(@(x)double(x),varargin,'UniformOutput',false);
            newobj = horzcat(c1{:});  % result is a double
        end
        function newobj = vertcat(varargin)
            c1 = cellfun(@(x)double(x),varargin,'UniformOutput',false);
            d = vertcat(c1{:});
            c2 = cellfun(@(x)x.absent,varargin,'UniformOutput',false);
            a = vertcat(c2{:});
            newobj = internal.stats.DoubleTableColumn(d,a);
        end
        function newobj = cat(dim,varargin)
            if isequal(dim,1)
                newobj = vertcat(varargin{:});
            else
                error(message('stats:internal:DoubleTableColumn:BadDim'));
            end
        end
    end
end
