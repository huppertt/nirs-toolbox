function newProps(d,h)
%NEWPROPS Create new properties in container object according to filter type.
%
%   Inputs:
%       d - Container object
%       h - Contained object


%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


% Get an array of structures with the property information
p = whichspecs(h);

% Create these properties in the container object 
createProps(d,h,p);

%--------------------------------------------------------------------
function createProps(d,h,p)
% Create these properties in the container object 
%
%
%   Inputs:
%       d - Container object
%       h - Contained object
%       p - array of properties to be replicated


% Delete existing dynamic properties
pold = get(d,'dynamicProps');
delete(pold);

% Clear listeners to existing dynamic properties
oldl = get(d,'dynamicPropsListener');
if isa(oldl,'handle.listener'), % Don't attempt delete if not a valid handle.listener
    delete(oldl);
end

% Do nothing if no new properties
if isempty(p),
    % Flush the property
    set(d,'dynamicProps',[]);
else
    lc = 0; % listener counter
    % Create dynamic properties that correspond to each property in p
    for n = 1:length(p),
        pnew(n) = schema.prop(d,p(n).name,p(n).datatype);
        set(d,p(n).name,p(n).defval);
        pnew(n).Description = p(n).descript;
        
        % Install a listener to the property when a callback is present
        if ~isempty(p(n).callback),
            lc = lc + 1;
            l(lc) = handle.listener(d,pnew(n),p(n).callback{1},p(n).callback{2});
            % Set callback target to object that generated the property
            cbt = findcbObj(h,p(n).descript);
            set(l(lc), 'callbacktarget', cbt); 
        end
    end
    
    % Store the handles to the dynamic properties
    set(d,'dynamicProps',pnew);
    
    % Store listener
    if lc > 0,
        set(d,'dynamicPropsListener',l);
    else
        % Flush the listeners
        set(d,'dynamicPropsListener',[]);
    end
end
