function [F, A] = getmask(this, fcns, rcf, specs)
%GETMASK   Get the mask.

%   Copyright 2008 The MathWorks, Inc.

specs = getspecs(this.CurrentSpecs);
     
PNF = specs.PeakNotchFrequencies(2:end);
units = feval(fcns.getunits);
LowerLimit = -60; 
%Draw a mask except when visualizing a zerophase response
if ~isequal(units,'zerophase')    
    if isequal(lower(this.Specification),'n,q')
        BW = (2/specs.FilterOrder)/specs.Q;
    else
        BW = specs.BW;
    end
    
    if isequal(lower(this.Specification),'l,bw,gbw,nsh')
        if this.ShelvingFilterOrder > 1
            LowerLimit = -120;
        end
    end
    
    Fq = [PNF-BW/2; PNF-BW/2; PNF-BW/2; PNF+BW/2; PNF+BW/2; PNF+BW/2];
    F = [0; BW/2; BW/2; BW/2; Fq(:)];
    A = [0 0 LowerLimit -inf -inf LowerLimit]';
    A = A(:,ones(1,floor(specs.FilterOrder/2)));
    A = [A(:); 0; 0];
    if isequal(lower(specs.CombType), 'peak')       
        if ~rem(specs.FilterOrder,2)
            F = [F(F<1); 1];            
        else
            F = [F; 1];
            A = [A; LowerLimit; -inf; -inf];
        end
    else
        if ~rem(specs.FilterOrder,2)
            F = [F(F<1); 1];                        
            A = [-inf; -inf; LowerLimit; A(1:length(F)-4); -inf];
        else
            F = [F; 1];
           A = [-inf; -inf; LowerLimit; A];
        end            
    end
        
    %scale frequency and amplitudes if necessary
    F = F*fcns.getfs()/2;
    
    if isequal(units,'linear')
        A = 10.^(A/20);
    end
    if isequal(units,'squared')
        A = 10.^(A/10);
    end       
else %zero phase response
    F=[];
    A=[];
end
end

% [EOF]
