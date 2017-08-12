function data = bemdata_spectral(mesh,frequency,wv)

% data = bemdata_spectral(mesh,frequency,wv)
%
% Calculates data (phase and amplitude) for a given
% spectral mesh at a given frequency (MHz).
% outputs phase and amplitude in structure data
% and mesh information in mesh
%
% mesh is the input mesh (variable or filename)
% frequency is the modulation frequency (MHz)
% wv is optional wavelength array



%% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

%% load mesh
if ischar(mesh)== 1
  mesh = load_mesh(mesh);
end

if nargin == 2
  wv_array = sort(mesh.wv);
elseif nargin == 3
  wv_array = sort(wv);
  for i = 1:length(wv_array)
    flag = find(mesh.wv == wv_array(i));
    if isempty(flag) == 1
      disp('ERROR: wv_array contains wavelengths not present in the extinction coefficients')
      data = [];
      return
    end
  end
end

nwv = length(wv_array);

%% loop through wavelengths
data.paa = []; data.wv=[];
for i = 1 : nwv
  disp(sprintf('Calculating data for: %g nm',(wv_array(i))))
  
  % calculate absorption and scattering coefficients from concetrations and
  % scattering parameters a and b
  [mesh.mua, mesh.mus, mesh.kappa] = calc_mua_mus(mesh,wv_array(i));
  
  % if sources are not fixed, move sources depending on mus
  if mesh.source.fixed == 0
    mus_eff = mesh.mus;
    [mesh]=move_source(mesh,mus_eff,3);
    clear mus_eff
  end
  
  [data_single_wv] = bemdata_stnd(mesh,frequency);
  data.paa = [data.paa, data_single_wv.paa];
  data.wv = [data.wv wv_array(i)];
  clear data_single_wv
end

data.link = mesh.link;
