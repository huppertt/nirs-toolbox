function [J,data,mesh] = jacobian(mesh,frequency,wv)

% [J,data,mesh]=jacobian(mesh,frequency,wv)
%
% Calculates Jacobian and data for a given mesh at given frequency 
% mesh can either be a workspace variable or filename.
% Automatically determines type.
% 
% mesh is the input mesh
% frequency is the modulation frequency (MHz)
% wv is optional wavelength array for spectral (e.g. [661 705 735])
% J is the Jacobian
% data is the field/boundary data



% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

% If not a workspace variable, load mesh
if ischar(mesh)== 1
    mesh = load_mesh(mesh);
end

% Select jacobian program based on mesh.type
if strcmp(mesh.type,'stnd') == 1        % stnd
  [J,data] = jacobian_stnd(mesh,frequency);
elseif strcmp(mesh.type,'stnd_bem')
    [J,data] = jacobian_stnd_bem(mesh,frequency); % bem stnd
elseif strcmp(mesh.type,'fluor') == 1   % fluor
  [data,mesh] = femdata_fl(mesh,frequency);
  datax.phi = data.phix;
  [J,data] = jacobian_fl(mesh,frequency,datax);
elseif strcmp(mesh.type,'fluor_bem') == 1   % fluor bem
  [data,mesh] = bemdata_fl(mesh,frequency);
  datax.phi = data.phix;
  [J,data] = jacobian_fl_bem(mesh,frequency,datax);
elseif strcmp(mesh.type,'spec') == 1    % spec
  if nargin == 2
      if frequency == 0
          [J,data] = jacobian_spectral_cw(mesh);
      else
        [J,data] = jacobian_spectral(mesh,frequency);
      end
  elseif nargin == 3
      if frequency == 0
          [J,data] = jacobian_spectral_cw(mesh,wv);
      else
        [J,data] = jacobian_spectral(mesh,frequency,wv);
      end
  end
elseif strcmp(mesh.type,'spec_bem') == 1    % spec bem
  if nargin == 2
      if frequency == 0
          [J,data] = jacobian_spectral_cw_bem(mesh);
      else
        [J,data] = jacobian_spectral_bem(mesh,frequency);
      end
  elseif nargin == 3
      if frequency == 0
          [J,data] = jacobian_spectral_cw_bem(mesh,wv);
      else
        [J,data] = jacobian_spectral_bem(mesh,frequency,wv);
      end
  end
end