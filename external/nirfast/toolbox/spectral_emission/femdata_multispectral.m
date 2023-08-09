function data = femdata_multispectral(mesh,frequency,drug,...
            wv_array_emiss,wv_excite,eta_tot_true,tau)
        
% data = femdata_multispectral(mesh,frequency,drug,...
%            wv_array_emiss,wv_excite,eta_tot_true,tau)  
%
% Forward solver - Calculates boundary data for a given mesh, 
% in both fluorescence and spectral cases for all wavelengths
%
% mesh is the input mesh (struct variable or filename)
% frequency is the modulation frequency (MHz)
% drug is the reemission spectrum (filename)
% wv_array_emiss is the wavelength array for emission
% wv_excite is the excitation wavelength
% eta_tot_true is the total quantum yield (scalar)
% tau is the lifetime (scalar)
        


%% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

%% load mesh
if ischar(mesh) == 1
    mesh = load_mesh(mesh);
end
if ~strcmp(mesh.type,'spec')
    errordlg('Mesh type must be spectral (''spec'')','NIRFAST Error');
    error('Mesh type must be spectral (''spec'')');
end

%% load emission spectrum
if ~exist(drug,'file')
    errordlg('Drug file not found','NIRFAST Error');
    error('Drug file not found');
end
emiss_spec_temp = load(drug);
emiss_spec.wv = emiss_spec_temp(:,1);
emiss_spec.etaspec = emiss_spec_temp(:,2);
clear emiss_spec_temp;

%% determine excitation wavelength parameters
[mesh.muax, mesh.musx, mesh.kappax, junk] = calc_mua_mus(mesh,wv_excite);
ex = mesh.excoef(find(mesh.wv == wv_excite),end);
[mesh.muaf] = mesh.conc(:,end)*ex;

data.paa = [];
data.wv = wv_array_emiss;

%% find emission field at each wavelength
for i = 1:numel(wv_array_emiss)
    disp(sprintf('Calculate fluorescence emission field at %g nm',wv_array_emiss(i))); 
    % determine mua mus values from concetrations and scatt params
    [mesh.muam, mesh.musm, mesh.kappam, junk] = calc_mua_mus(mesh,wv_array_emiss(i));
    
    emiss_wv_min = min(abs(emiss_spec.wv - wv_array_emiss(i))); 
    emiss_diff = abs(emiss_spec.wv - wv_array_emiss(i));
    emiss_wv_index = find(emiss_diff == emiss_wv_min(1));
    mesh.eta = repmat(emiss_spec.etaspec(emiss_wv_index(1))*eta_tot_true,...
        size(mesh.muam,1),1);
    mesh.tau = zeros(size(mesh.muam,1),1);
    mesh.tau(:) = tau;

    data_tmp = femdata_fl(mesh,frequency);
    data.paa = [data.paa data_tmp.paafl];
end

%% put link info into data
data.link = mesh.link;

