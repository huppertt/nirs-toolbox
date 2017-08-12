function [data,mesh] = femdata(mesh,frequency,wv)

% [data,mesh]=femdata(mesh,frequency,wv)
%
% Calculates forward model data for a given mesh at given frequency 
%
% mesh can either be a workspace variable or filename
%   (automatically determines type)
% frequency is the modulation frequency (MHz)
% wv is an optional wavelength array for spectral (e.g. [661 705 735])


% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

% If not a workspace variable, load mesh
if ischar(mesh)== 1
  mesh = load_mesh(mesh);
end

% Select femdata program based on mesh.type
if strcmp(mesh.type,'stnd') == 1        % Use stnd NIRFAST, non-spectral mesh
  data = femdata_stnd(mesh,frequency);
elseif strcmp(mesh.type,'fluor') == 1   % Use stnd fluorfast, non-spectral mesh
  data = femdata_fl(mesh,frequency);
elseif strcmp(mesh.type,'spec') == 1    % Use spectral mesh
  if nargin == 2
    data = femdata_spectral(mesh,frequency);
  elseif nargin == 3
    data = femdata_spectral(mesh,frequency,wv);
  end
elseif strcmp(mesh.type,'stnd_spn') == 1        % stnd spn (assume n=5 if this was called)
  data = femdata_sp5(mesh,frequency);
elseif strcmp(mesh.type,'stnd_bem') == 1        % stnd bem
  data = bemdata_stnd(mesh,frequency);
elseif strcmp(mesh.type,'fluor_bem') == 1   % fluor bem
  data = bemdata_fl(mesh,frequency);
elseif strcmp(mesh.type,'spec_bem') == 1    % spec bem
  if nargin == 2
    data = bemdata_spectral(mesh,frequency);
  elseif nargin == 3
    data = bemdata_spectral(mesh,frequency,wv);
  end
else
    errordlg('The mesh type was not found; you may need to call a specific forward solver','NIRFAST Error');
    error('The mesh type was not found; you may need to call a specific forward solver');
end
