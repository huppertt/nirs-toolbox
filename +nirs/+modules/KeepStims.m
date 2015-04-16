classdef KeepStims < nirs.modules.AbstractModule
    properties
        listOfStims = {};
    end
    
    properties( Access = private )
       discard
    end
    
    methods
        function obj = KeepStims( prevJob )
           obj.name     = 'Keep Just These Stims';
           obj.discard  = nirs.modules.DiscardStims( prevJob );
           
           obj.discard.keepOrDiscard  = 'keep';
        end
        
        function data = runThis( obj, data )
            obj.discard.listOfStims = obj.listOfStims;
            data = obj.discard.run( data );
        end
    end
    
end

