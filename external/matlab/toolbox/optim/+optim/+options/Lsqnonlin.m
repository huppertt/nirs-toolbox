classdef (Sealed) Lsqnonlin < optim.options.Lsqncommon
%

%Lsqnonlin Options for LSQNONLIN
%
%   The OPTIM.OPTIONS.LSQNONLIN class allows the user to create a set of
%   options for the LSQNONLIN solver. For a list of options that can be set,
%   see the documentation for LSQNONLIN.
%
%   OPTS = OPTIM.OPTIONS.LSQNONLIN creates a set of options for LSQNONLIN
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.LSQNONLIN(PARAM, VAL, ...) creates a set of options
%   for LSQNONLIN with the named parameters altered with the specified
%   values.
%
%   OPTS = OPTIM.OPTIONS.LSQNONLIN(OLDOPTS, PARAM, VAL, ...) creates a copy
%   of OLDOPTS with the named parameters altered with the specified values.
%
%   See also OPTIM.OPTIONS.MULTIALGORITHM, OPTIM.OPTIONS.SOLVEROPTIONS

%   Copyright 2012-2015 The MathWorks, Inc.    
       
    properties (Hidden)
%SOLVERNAME Name of the solver that the options are intended for
%         
        SolverName = 'lsqnonlin';
    end       

    properties (Hidden, SetAccess = private, GetAccess = public)
        
        % New version property added in third version
        LsqnonlinVersion
    end
    
    methods (Hidden)
        
        function obj = Lsqnonlin(varargin)
%Lsqnonlin Options for LSQNONLIN
%
%   The OPTIM.OPTIONS.LSQNONLIN class allows the user to create a set of
%   options for the LSQNONLIN solver. For a list of options that can be set,
%   see the documentation for LSQNONLIN.
%
%   OPTS = OPTIM.OPTIONS.LSQNONLIN creates a set of options for LSQNONLIN
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.LSQNONLIN(PARAM, VAL, ...) creates a set of options
%   for LSQNONLIN with the named parameters altered with the specified
%   values.
%
%   OPTS = OPTIM.OPTIONS.LSQNONLIN(OLDOPTS, PARAM, VAL, ...) creates a copy
%   of OLDOPTS with the named parameters altered with the specified values.
%
%   See also OPTIM.OPTIONS.MULTIALGORITHM, OPTIM.OPTIONS.SOLVEROPTIONS

            % Call the superclass constructor
            obj = obj@optim.options.Lsqncommon(varargin{:});
               
            % Record the class version; Update property 'LsqnonlinVersion'
            % instead of superclass property 'Version'.
            obj.Version = 2;
            obj.LsqnonlinVersion = 3;    
        end        
        
    end

    methods (Static)
        
        function obj = loadobj(obj)
    
            % Objects saved in R2013a will come in as structures. 
            if isstruct(obj) && obj.Version == 1

                % Save the existing structure
                s = obj;
                
                % Create a new object
                obj = optim.options.Lsqnonlin;
                
                % Call the superclass method to upgrade the object
                obj = upgradeFrom13a(obj, s); 
                
                % The SolverVersion property was not present in 13a. We
                % clear it here and the remainer of loadobj will set it
                % correctly.
                obj.LsqnonlinVersion = [];
                
            end
            
            % Upgrading to 13b
            % Update Levenberg-Marquardt and finite difference options
            if obj.Version < 2
                obj = upgradeLMandFinDiffTo13b(obj);
            end
                
            % Set the superclass version number
            obj.Version = 2;

            % Upgrading to 15b            
            % Introduce LsqnonlinVersion field
            if isempty(obj.LsqnonlinVersion)
                % Use 'LsqnonlinVersion' property instead of 'Version'
                % property because 'Version' is for the superclass and
                % 'LsqnonlinVersion' is for this (derived) class. However,
                % 'LsqnonlinVersion' was added in the second version, we
                % check only for the second version and add this property.
                % For all other version, check only the 'LsqnonlinVersion'
                % property.
                obj.LsqnonlinVersion = 2; % update object
            end
               
            % Set the version number
            obj.LsqnonlinVersion = 3;            
        end
        
    end
            
    
end
    