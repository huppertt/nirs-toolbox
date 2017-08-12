classdef titleddataset < dataset
        
%   Copyright 2012-2014 The MathWorks, Inc.    
    
    properties (Access=private)
        Title = '';
    end
    
    methods (Hidden=true)
        function ds = titleddataset(varargin)
            if isa(varargin{1},'table')
                varargin{1} = table2dataset(varargin{1});
            end
            copying = nargin>=1 && isa(varargin{1},'dataset');
            if copying
                args = {};
            else
                args = varargin;
            end
            ds@dataset(args{:});
            if copying
                ds = classreg.regr.lmeutils.datasetcopier(ds,varargin{1});
                if nargin>=2
                    ds.Title = varargin{2};
                end
            end
        end
    end
    
    methods
        function ds = settitle(ds,title)
            ds.Title = title;
        end
        function disp(ds)
            if ~isempty(ds.Title)
                if feature('hotlinks')
                    fprintf('\n    <strong>%s</strong>\n\n',ds.Title)
                else
                    fprintf('\n    %s\n\n',upper(ds.Title))
                end
            end
          disp@dataset(ds)
        end
    end
    
end

