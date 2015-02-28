function mesh = set_mesh_type(mesh,type)

% mesh = set_mesh_type(mesh,type)
%
% changes mesh type and optical properties
%
% mesh is the variable or mesh location
% type is the mesh type


%% make sure it's not a conversion between fem and bem
if isfield(mesh,'type')
    bem1 = 0;
    if strcmp(mesh.type,'stnd_bem') || strcmp(mesh.type,'fluor_bem') ...
            || strcmp(mesh.type,'spec_bem')
        bem1 = 1;
    end
    bem2 = 0;
    if strcmp(type,'stnd_bem') || strcmp(type,'fluor_bem') ...
            || strcmp(type,'spec_bem')
        bem2 = 1;
    end
    
    if bem1 + bem2 == 1
        errordlg('Cannot convert between FEM and BEM','NIRFAST Error');
        error('Cannot convert between FEM and BEM');
    end
end

%% assign properties based on type
if strcmp(type,'stnd')
    mesh.mua = ones(length(mesh.nodes),1).*0.01;
    mesh.mus = ones(length(mesh.nodes),1);
    mesh.kappa = 1./(3.*(mesh.mua+mesh.mus));
    mesh.ri = ones(length(mesh.nodes),1).*1.33;
elseif strcmp(type,'fluor')
    mesh.muax = ones(length(mesh.nodes),1).*0.009;
    mesh.musx = ones(length(mesh.nodes),1).*1.314;
    mesh.kappax = 1./(3.*(mesh.muax+mesh.musx));
    mesh.muam = ones(length(mesh.nodes),1).*0.006;
    mesh.musm = ones(length(mesh.nodes),1).*1.273;
    mesh.kappam = 1./(3.*(mesh.muam+mesh.musm));
    mesh.muaf = ones(length(mesh.nodes),1).*0.002;
    mesh.eta = ones(length(mesh.nodes),1).*0.1;
    mesh.tau = ones(length(mesh.nodes),1).*0;
    mesh.ri = ones(length(mesh.nodes),1).*1.33;
elseif strcmp(type,'spec')
    mesh.sa = ones(length(mesh.nodes),1);
    mesh.sp = ones(length(mesh.nodes),1);
    c1 = ones(length(mesh.nodes),1).*0.01;
    c2 = ones(length(mesh.nodes),1).*0.01;
    c3 = ones(length(mesh.nodes),1).*0.4;
    mesh.conc = [c1 c2 c3];
    mesh.chromscattlist = [{'HbO'};{'deoxyHb'};{'Water'};{'S-Amplitude'};{'S-Power'}];
    mesh.wv = [661;735;761;785;808;826;849];
    mesh.excoef = [0.0741    0.8500    0.0015;
                    0.0989    0.2400    0.0038;
                    0.1185    0.3292    0.0043;
                    0.1500    0.2056    0.0038;
                    0.1741    0.1611    0.0033;
                    0.2278    0.1611    0.0040;
                    0.2370    0.1556    0.0058];
    if isfield(mesh,'link')
        mesh.link(:,end+1:end+6) = 1;
    end
elseif strcmp(type,'stnd_spn')
    mesh.mua = ones(length(mesh.nodes),1).*0.01;
    mesh.mus = ones(length(mesh.nodes),1);
    mesh.kappa = 1./(3.*(mesh.mua+mesh.mus));
    mesh.ri = ones(length(mesh.nodes),1).*1.33;
    mesh.g = ones(length(mesh.nodes),1).*0.9;
elseif strcmp(type,'stnd_bem')
    nregions = size(unique(mesh.region),1)-1;
    mesh.mua = ones(nregions,1).*0.01;
    mesh.mus = ones(nregions,1);
    mesh.kappa = 1./(3.*(mesh.mua+mesh.mus));
    mesh.ri = ones(nregions,1).*1.33;
elseif strcmp(type,'fluor_bem')
    nregions = size(unique(mesh.region),1)-1;
    mesh.muax = ones(nregions,1).*0.009;
    mesh.musx = ones(nregions,1).*1.314;
    mesh.kappax = 1./(3.*(mesh.muax+mesh.musx));
    mesh.muam = ones(nregions,1).*0.006;
    mesh.musm = ones(nregions,1).*1.273;
    mesh.kappam = 1./(3.*(mesh.muam+mesh.musm));
    mesh.muaf = ones(nregions,1).*0.002;
    mesh.eta = ones(nregions,1).*0.1;
    mesh.tau = ones(nregions,1).*0;
    mesh.ri = ones(nregions,1).*1.33;
elseif strcmp(type,'spec_bem')
    nregions = size(unique(mesh.region),1)-1;
    mesh.sa = ones(nregions,1);
    mesh.sp = ones(nregions,1);
    c1 = ones(nregions,1).*0.01;
    c2 = ones(nregions,1).*0.01;
    c3 = ones(nregions,1).*0.4;
    mesh.conc = [c1 c2 c3];
    mesh.chromscattlist = [{'HbO'};{'deoxyHb'};{'Water'};{'S-Amplitude'};{'S-Power'}];
    mesh.wv = [661;735;761;785;808;826;849];
    mesh.excoef = [0.0741    0.8500    0.0015;
                    0.0989    0.2400    0.0038;
                    0.1185    0.3292    0.0043;
                    0.1500    0.2056    0.0038;
                    0.1741    0.1611    0.0033;
                    0.2278    0.1611    0.0040;
                    0.2370    0.1556    0.0058];
else
    errordlg('Invalid mesh type','NIRFAST Error');
    error('Invalid mesh type: %s', type);
end

mesh.type = type;
