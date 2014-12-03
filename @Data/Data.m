classdef Data
    %DATA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        time;
        data;
        probe;
        modFreq;
        description;
    end
    
    methods
        %% Constructor
        function obj = Data( varargin )
            
            if nargin > 0
                obj.time = varargin{1};
            end
            
            if nargin > 1
                obj.data = varargin{2};
            end
            
            if nargin > 2
                obj.probe = varargin{3};
            end
            
            if nargin > 3
                obj.modFreq = varargin{4};
            end
            
            if nargin > 4
                obj.description = varargin{5};
            end
            
            if nargin > 5
                error('Too many input arguments.')
            end
        end
        
        %% Set/Get
        function obj = set.time( obj, newTime )
            if iscolumn( newTime )
                obj.time = newTime;
            elseif ~iscolumn( newTime ) && isrow( newTime )
                obj.time = newTime.';
            else
                error( 'Time must be a column vector.')
            end     
        end

        function obj = set.data( obj, newData )
            if ismatrix( newData ) && isnumeric( newData );
                obj.data = newData;
            end
        end
        
        function obj = set.modFreq( obj, newFreq )
            if newFreq >= 0 && numel( newFreq ) == 1
                obj.modFreq = newFreq;
            else
                error( 'Modulation frequency should not be negative.')
            end
        end
        
        function obj = set.probe( obj, newProbe )
            if strcmpi( class(newProbe),'nirs.Probe' )
                obj.probe = newProbe;
            else
                error('Probe must be of class "Probe".')
            end
        end
        
        function obj = set.description( obj, newDescription )
            if ischar( newDescription ) || isempty( newDescription )
                obj.description = newDescription;
            else
                error( 'Description must be a string.' )
            end
        end
    end
end

