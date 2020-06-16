function mesh = set_mesh(mesh, region, values)

% mesh = set_mesh(mesh, region, values)
%
% assigns the values to the given region in the mesh
%
% mesh is the mesh filename or variable
% region is the region number
% values is a structure that contains the desired values
%
% for a STANDARD mesh:
%  values.mua - absorption coefficient
%  values.mus - scatter coefficient
%  values.ri - refractive index
%
% for a STANDARD BEM mesh:
%  values.mua - absorption coefficient
%  values.mus - scatter coefficient
%  values.ri - refractive index
%
% for a STANDARD SPN mesh:
%  values.mua - absorption coefficient
%  values.mus - scatter coefficient
%  values.g - anisotropic number
%  values.ri - refractive index
% 
% for a FLUORESCENCE or FLUORESCENCE BEM mesh:
%  values.muax - absorption coefficient at excitation
%  values.musx - scatter coefficient at excitation
%  values.muam - absorption coefficient at emission
%  values.musm - scatter coefficient at emission
%  values.muaf - absorption of fluorophore
%  values.eta - quantum yield of fluorophore
%  values.tau - lifetime of fluorophore
%  values.ri - refractive index
% 
% for a SPECTRAL or SPECTRAL BEM mesh:
%  values.sa - S-amplitude
%  values.sp - S-power
%  values.HbO, values.deoxyHb, etc. for each chromophore
%  values.ri - refractive index



% If not a workspace variable, load mesh
if ischar(mesh)== 1
  mesh = load_mesh(mesh);
end

if ~strcmp(mesh.type,'stnd_bem')
    ind = find(mesh.region==region);
end

% STANDARD
if strcmp(mesh.type,'stnd')
    if isfield(values,'mua')
        mesh.mua(ind) = values.mua;
        mesh.kappa(ind) = 1./(3.*(mesh.mua(ind)+mesh.mus(ind)));
    end
    if isfield(values,'mus')
        mesh.mus(ind) = values.mus;
        mesh.kappa(ind) = 1./(3.*(mesh.mua(ind)+mesh.mus(ind)));
    end
    if isfield(values,'ri')
        mesh.ri(ind) = values.ri;
        mesh.c(ind)=(3e11/values.ri);
    end

% STANDARD BEM
elseif strcmp(mesh.type,'stnd_bem')
    if isfield(values,'mua')
        mesh.mua(region) = values.mua;
        mesh.kappa(region) = 1./(3.*(mesh.mua(region)+mesh.mus(region)));
    end
    if isfield(values,'mus')
        mesh.mus(region) = values.mus;
        mesh.kappa(region) = 1./(3.*(mesh.mua(region)+mesh.mus(region)));
    end
    if isfield(values,'ri')
        mesh.ri(region) = values.ri;
        mesh.c(region)=(3e11/values.ri);
    end    
    
% STANDARD SPN
elseif strcmp(mesh.type,'stnd_spn')
    if isfield(values,'mua')
        mesh.mua(ind) = values.mua;
        mesh.kappa(ind) = 1./(3.*(mesh.mua(ind)+mesh.mus(ind)));
    end
    if isfield(values,'mus')
        mesh.mus(ind) = values.mus;
        mesh.kappa(ind) = 1./(3.*(mesh.mua(ind)+mesh.mus(ind)));
    end
    if isfield(values,'g')
        mesh.g(ind) = values.g;
    end
    if isfield(values,'ri')
        mesh.ri(ind) = values.ri;
        mesh.c(ind)=(3e11/values.ri);
    end
    
% FLUORESCENCE BEM
elseif strcmp(mesh.type,'fluor_bem')
    if isfield(values,'muax')
        mesh.muax(region) = values.muax;
        mesh.kappax(region) = 1./(3.*(mesh.muax(region)+mesh.musx(region)));
    end
    if isfield(values,'musx')
        mesh.musx(region) = values.musx;
        mesh.kappax(region) = 1./(3.*(mesh.muax(region)+mesh.musx(region)));
    end
    if isfield(values,'muam')
        mesh.muam(region) = values.muam;
        mesh.kappam(region) = 1./(3.*(mesh.muam(region)+mesh.musm(region)));
    end
    if isfield(values,'musm')
        mesh.musm(region) = values.musm;
        mesh.kappam(region) = 1./(3.*(mesh.muam(region)+mesh.musm(region)));
    end
    if isfield(values,'muaf')
        mesh.muaf(region) = values.muaf;
    end
    if isfield(values,'eta')
        mesh.eta(region) = values.eta;
    end
    if isfield(values,'tau')
        mesh.tau(region) = values.tau;
    end
    if isfield(values,'ri')
        mesh.ri(region) = values.ri;
        mesh.c(region)=(3e11/values.ri);
    end
    
% FLUORESCENCE
elseif strcmp(mesh.type,'fluor')
    if isfield(values,'muax')
        mesh.muax(ind) = values.muax;
        mesh.kappax(ind) = 1./(3.*(mesh.muax(ind)+mesh.musx(ind)));
    end
    if isfield(values,'musx')
        mesh.musx(ind) = values.musx;
        mesh.kappax(ind) = 1./(3.*(mesh.muax(ind)+mesh.musx(ind)));
    end
    if isfield(values,'muam')
        mesh.muam(ind) = values.muam;
        mesh.kappam(ind) = 1./(3.*(mesh.muam(ind)+mesh.musm(ind)));
    end
    if isfield(values,'musm')
        mesh.musm(ind) = values.musm;
        mesh.kappam(ind) = 1./(3.*(mesh.muam(ind)+mesh.musm(ind)));
    end
    if isfield(values,'muaf')
        mesh.muaf(ind) = values.muaf;
    end
    if isfield(values,'eta')
        mesh.eta(ind) = values.eta;
    end
    if isfield(values,'tau')
        mesh.tau(ind) = values.tau;
    end
    if isfield(values,'ri')
        mesh.ri(ind) = values.ri;
        mesh.c(ind)=(3e11/values.ri);
    end

% SPECTRAL BEM
elseif strcmp(mesh.type,'spec_bem')
    if isfield(values,'sa')
        mesh.sa(region) = values.sa;
    end
    if isfield(values,'sp')
        mesh.sp(region) = values.sp;
    end
    if isfield(values,'ri')
        mesh.ri(region) = values.ri;
        mesh.c(region)=(3e11/values.ri);
    end
    for i=1:1:numel(mesh.chromscattlist)-2
        if isfield(values,mesh.chromscattlist(i))
            fld=mesh.chromscattlist(i);
            comm = strcat('mesh.conc(region,i)=values.',fld,';');
            eval(comm{1});
        end
    end
    
% SPECTRAL
elseif strcmp(mesh.type,'spec')
    if isfield(values,'sa')
        mesh.sa(ind) = values.sa;
    end
    if isfield(values,'sp')
        mesh.sp(ind) = values.sp;
    end
    if isfield(values,'ri')
        mesh.ri(ind) = values.ri;
        mesh.c(ind)=(3e11/values.ri);
    end
    for i=1:1:numel(mesh.chromscattlist)-2
        if isfield(values,mesh.chromscattlist(i))
            fld=mesh.chromscattlist(i);
            comm = strcat('mesh.conc(ind,i)=values.',fld,';');
            eval(comm{1});
        end
    end
end