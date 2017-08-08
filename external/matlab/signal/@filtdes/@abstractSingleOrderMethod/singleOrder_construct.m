function singleOrder_construct(d,varargin)
%SINGLEORDER_CONSTRUCT  'Real' constructor for the single order design method object.


%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.



% Create a dynamic property to hold the filter order
% Use a dynamic property so that its visibility can
% be turned off in the min order class that inherits
% from this one
p = schema.prop(d,'order','spt_uint32');
set(d,'order',8);

% Call super's constructor, do this after adding the order prop
designMethodwFs_construct(d,varargin{:});

