classdef KeepStims < nirs.modules.AbstractModule
    properties
        list = {};
    end
    
    properties( Access = private )
       discard
    end
    
    methods
        function obj = KeepStims( prevJob )
           obj.name     = 'Keep Just These Stims';
           obj.discard  = nirs.modules.DiscardStims( prevJob );
           obj.flag     = 'keep';
        end
        
        function data = runThis( obj, data )
            data = obj.discard.run( data );
        end
    end
    
end

