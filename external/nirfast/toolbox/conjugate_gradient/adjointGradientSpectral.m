function [ z, data, z_amp, z_phase ] = adjointGradientSpectral( mesh, measurements,frequency,regParam )
%ADJOINTGRADIENTSPECTRAL Calculates gradient for spectral meshes using
%   adjoint method.
%% If no regularisation parameters have been provided, assume none.
if nargin < 4
    regParam.lamba = 0;
end
if ~isfield(mesh,'cgFormatLink')
    mesh = convertLinkFormat(mesh);
end
%% Preallocate gradient and data vectors.
nnodes = size(mesh.nodes,1);
nwv = length(mesh.wv);
npara = size(mesh.conc,2);
z = zeros(nnodes*(npara+2),1);
z_amp = zeros(nnodes*(npara+2),1);
z_phase = zeros(nnodes*(npara+2),1);
data.paa = zeros(size(mesh.cgFormatLink,1),nwv*2);
%% For each wavelength, solve standard problem and use solution to
%% calculate gradients.
mesh.type = 'stnd';
for ii = 1:nwv
    tempMesh = mesh;
    %Calculate mesh optical parameters and move sources.
    [tempMesh.mua, tempMesh.mus, tempMesh.kappa] = calc_mua_mus(tempMesh,tempMesh.wv(ii));
    
    % if sources are not fixed, move sources based on mus.
    if tempMesh.source.fixed == 0
        mus_eff = tempMesh.mus;
        [tempMesh]=move_source(tempMesh,mus_eff,3);
        clear mus_eff
    end
    templink = zeros(size(mesh.cgFormatLink,1),3);
    templink(:,1:2) = mesh.cgFormatLink(:,1:2);
    templink(:,3) = mesh.cgFormatLink(:,2+ii);
    tempMesh.link = templink;
    tempMesh.cgFormatLink(:,3) = templink(:,3);
    %% Calculate gradient for optical parameters and current wavelength.
    tempMeas = struct;
    tempMeas.paa = measurements.paa(:,(2*ii-1):(2*ii));
    if iscell(tempMesh.R)
        tempMesh.R = tempMesh.R{ii};
    end
    tempRegParam = regParam;
    tempRegParam.P = regParam.P{ii};
    [z_i data_i z_i_amp z_i_phase] = adjointGradientStnd(tempMesh,tempMeas,frequency,tempRegParam);
    z_i_kappa = z_i(1:nnodes);
    z_i_mua = z_i((nnodes+1):end);
    z_i_kappa_amp = z_i_amp(1:nnodes);
    z_i_kappa_phase = z_i_phase(1:nnodes);
    z_i_mua_amp = z_i_amp((nnodes+1):end);
    z_i_mua_phase = z_i_phase((nnodes+1):end);
    
    %Record forward data for current wavelength.
    data.paa(:,(2*ii-1):(2*ii)) = data_i.paa;
    %Calculate concentration gradients (extinction
    %coefficient*d_phi/d_mua).
    for jj = 1:npara
        startInd = (jj-1)*nnodes+1;
        endInd = jj*nnodes;
        z(startInd:endInd) = ...
            z(startInd:endInd) + z_i_mua*tempMesh.excoef(ii,jj);
        z_amp(startInd:endInd) = ...
            z_amp(startInd:endInd) + z_i_mua_amp*tempMesh.excoef(ii,jj);
        z_phase(startInd:endInd) = ...
            z_phase(startInd:endInd) + z_i_mua_phase*tempMesh.excoef(ii,jj);
    end
    
    %Calculate scattering parameter gradients.
    dKappadMus = -3*(tempMesh.kappa.^2);
    dPhidMus = z_i_kappa.*dKappadMus;
    dPhidMus_amp = z_i_kappa_amp.*dKappadMus;
    dPhidMus_phase = z_i_kappa_phase.*dKappadMus;
    saStart = npara*nnodes+1;
    saEnd = (npara+1)*nnodes;
    
    %sa gradient is d_phi/d_mus * d_mus/d_sa.
    z(saStart:saEnd) = z(saStart:saEnd) + dPhidMus.*((tempMesh.wv(ii)/1000).^(-tempMesh.sp));
    z_phase(saStart:saEnd) = z_phase(saStart:saEnd) + dPhidMus_phase.*((tempMesh.wv(ii)/1000).^(-tempMesh.sp));
    z_amp(saStart:saEnd) = z_amp(saStart:saEnd) + dPhidMus_amp.*((tempMesh.wv(ii)/1000).^(-tempMesh.sp));
    
    spStart = (npara+1)*nnodes+1;
    spEnd = (npara+2)*nnodes;
    
    %sp gradient is d_phi/d_mus * d_mus/d_sp.
    z(spStart:spEnd) = z(spStart:spEnd) + dPhidMus.*(-tempMesh.sa.*((tempMesh.wv(ii)/1000).^(-tempMesh.sp))).*log(tempMesh.wv(ii)/1000); 
    z_phase(spStart:spEnd) = z_phase(spStart:spEnd) + dPhidMus_phase.*(-tempMesh.sa.*((tempMesh.wv(ii)/1000).^(-tempMesh.sp))).*log(tempMesh.wv(ii)/1000); 
    z_amp(spStart:spEnd) = z_amp(spStart:spEnd) + dPhidMus_amp.*(-tempMesh.sa.*((tempMesh.wv(ii)/1000).^(-tempMesh.sp))).*log(tempMesh.wv(ii)/1000); 
end

end

