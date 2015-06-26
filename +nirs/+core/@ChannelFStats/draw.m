function draw( obj, frange, idx )
    F       = obj.F;
    fcrit   = obj.Fcrit;

    if nargin < 2
        frange = [0 ceil(max(abs(F(:))))];
    end

    % loop through var names
    h = []; % handles
    types = obj.probe.link.type;

    if any(isnumeric(types))
        types = cellfun(@(x){num2str(x)}, num2cell(types));
    end

    utypes = unique(types, 'stable');

    for iName = 1:length( obj.names )

        for iType = 1:length(utypes)
            lst = strcmp( types, utypes(iType) );

            h(end+1) = figure;
            f = F(iName, lst);
            obj.probe.draw( f, frange, fcrit(iName) );
            title([utypes(iType) ' : ' obj.names{iName}], 'Interpreter','none')
        end
    end
end
