function s = obj2struct(this)
%OBJ2STRUCT <short description>

%   Copyright 2010 The MathWorks, Inc.

s = get(this);
if ~isempty(this.inport)
    for I = 1:length(this.inport)
        X(I) = obj2struct(this.inport(I));
    end
    s.inport = X;
end

clear X;

if ~isempty(this.outport)
    for I = 1:length(this.outport)
        X(I) = obj2struct(this.outport(I));
    end
    s.outport = X;
end

