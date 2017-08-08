classdef BlankGenerator < classreg.learning.generator.Generator

%   Copyright 2010-2011 The MathWorks, Inc.


    methods
        function this = BlankGenerator(X,Y,W,fitData)
            this = this@classreg.learning.generator.Generator(X,Y,W,fitData);
        end
    end

    methods
        function [this,X,Y,W,fitData,optArgs] = generate(this)
            X = this.X;
            Y = this.Y;
            W = this.W;
            fitData = this.FitData;
            optArgs = {};
        end
        
        function this = update(this,X,Y,W,fitData)
            this.X = X;
            this.Y = Y;
            this.W = W;
            this.FitData = fitData;
        end
     end

end
