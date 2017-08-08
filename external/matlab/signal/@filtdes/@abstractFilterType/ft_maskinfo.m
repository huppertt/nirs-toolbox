function cmd = ft_maskinfo(hObj, d)
%FT_MASKINFO Returns the mask information.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

specobjs = get(hObj, 'SpecObjs');

endx = [];

cmd = base_maskinfo(hObj, d);

cmd.bands = maskinfo(specobjs(1), d);
for indx = 2:length(specobjs),
    newcmd = maskinfo(specobjs(indx), d);
    
    % If the first specobj returned no commands, ignore it and just use the
    % 2nd specobj commands.
    if isempty(cmd.bands),
        cmd.bands = newcmd;
    elseif ~isempty(newcmd),
        
        % If the 2nd specobj returned commands integrate them
        for jndx = 1:length(cmd.bands),
            
            % If the index is empty ignore that band
            if isempty(newcmd{jndx}) || isempty(cmd.bands{jndx}),
                endx = [endx jndx];
            else
                cmd.bands{jndx} = setstructfields(cmd.bands{jndx}, newcmd{jndx});
            end
        end
    end
end

% Throw out all the empty indexes
cmd.bands(endx) = [];

% [EOF]
