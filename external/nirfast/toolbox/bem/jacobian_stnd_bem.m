function [J,data] = jacobian_stnd_bem(mesh,frequency)

% [J,data] = jacobian_stnd_bem(mesh,frequency)
%
% Calculates bem jacobian using perturbation method
%
% mesh is the filename or variable mesh
% frequency is the modulation frequency
% J is the resulting Jacobian
% data is the boundary data


%% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

%% load mesh
if ischar(mesh)== 1
  mesh = load_mesh(mesh);
end

%% calulate data
data = bemdata_stnd(mesh,frequency);

data.paa(:,1) = log(data.paa(:,1));
if frequency ~=0
    data.paa(:,2) = data.paa(:,2).*pi/180.0;
    data.paa(data.paa(:,2)<0,2) = data.paa(data.paa(:,2)<0,2) + (2*pi);
    data.paa(data.paa(:,2)>(2*pi),2) = data.paa(data.paa(:,2)>(2*pi),2) - (2*pi);
end
data.link = mesh.link;

%% jacobian
nregions = length(mesh.mua);
if frequency == 0
    J = zeros(length(data.paa), nregions);
else
    J = zeros(2*length(data.paa), 2*nregions);
end
mesh2 = mesh;

for reg = 1:nregions
    % mua
    mesh2.mua = mesh.mua;
    delta_mua = 0.01*mesh.mua(reg);
    mesh2.mua(reg) = mesh.mua(reg) + delta_mua;
    data2 = bemdata_stnd(mesh2,frequency);
    data2.paa(:,1) = log(data2.paa(:,1));
    if frequency ~= 0
        data2.paa(:,2) = data2.paa(:,2).*pi/180.0;
        data2.paa(data2.paa(:,2)<0,2) = data2.paa(data2.paa(:,2)<0,2) + (2*pi);
        data2.paa(data2.paa(:,2)>(2*pi),2) = data2.paa(data2.paa(:,2)>(2*pi),2) - (2*pi);
    end
    data_diff = data2.paa - data.paa;
    if frequency ~= 0
        data_diff = reshape(data_diff',2*length(data_diff),1);
    else
        data_diff = data_diff(:,1);
    end
    
    J(:,reg) = (data_diff)./(delta_mua);
end

if frequency ~= 0
    mesh2=mesh;
    for reg = 1:nregions
        % kappa
        mesh2.kappa = mesh.kappa;
        delta_kappa = 0.01*mesh.kappa(reg);
        mesh2.kappa(reg) = mesh.kappa(reg) + delta_kappa;
        data2 = bemdata_stnd(mesh2,frequency);
        data2.paa(:,1) = log(data2.paa(:,1));
        data2.paa(:,2) = data2.paa(:,2).*pi/180.0;
        data2.paa(data2.paa(:,2)<0,2) = data2.paa(data2.paa(:,2)<0,2) + (2*pi);
        data2.paa(data2.paa(:,2)>(2*pi),2) = data2.paa(data2.paa(:,2)>(2*pi),2) - (2*pi);
        data_diff = data2.paa - data.paa;
        data_diff = reshape(data_diff',2*length(data_diff),1);

        J(:,reg+nregions) = (data_diff)./(delta_kappa);
    end
end
