function newarith = set_privarith(this, newarith)
%SET_PRIVARITH   SetFunction for the 'privArithmetic' property.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

if isempty(newarith)
    newarith = defaultarithmetic(this);
end

constr = filtquant_plugins(this, newarith);

hfq = get(this, 'privfq');

if ~isempty(hfq)

    % Try to find object
    hfq = find(hfq,'-class',constr);
end

% Construct and store a new object if not found
if isempty(hfq),
    hfq = feval(constr);
    setdefaultcoeffwl(hfq,this);
    this.privfq = [this.privfq; hfq]; % Note that this.privfq(end+1) = hfq; doesn't work!
end

% Set the current filterquantizer
this.filterquantizer = hfq;

% [EOF]
