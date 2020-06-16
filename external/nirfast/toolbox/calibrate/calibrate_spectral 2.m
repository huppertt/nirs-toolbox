function [data,mesh,mua_big,mus_big] = calibrate_spectral(homog_data,...
                   anom_data,...
                   mesh_homog,...
                   mesh_anom,...
        		   frequency,...
        		   iteration,...
                   nographs,...
                   wv_array)

% [data,mesh,mua_big,mus_big] = calibrate_spectral(homog_data,...
%                   anom_data,...
%                   mesh_homog,...
%                   mesh_anom,...
%        		   frequency,...
%        		   iteration,...
%                   nographs,...
%                   wv_array)
%
% Calibrates spectral data, using homogeneous data measured on phantom
% and the anomaly data. This function requires the MATLAB optimization
% toolbox in order to make an initial guess.
%
% homog_data is the homogeneous data (filename or variable)
% anom_data is the anomaly data (filename or variable)
% mesh_homog is the homogeneous mesh (filename or variable)
% mesh_anom is the anomaly mesh (filename or variable)
% frequency is the modulation frequency (MHz)
% iteration is the number of iterations for fitting
% nographs is a flag for displaying the graphs
% wv_array is the desired wavelength array
% data is the calibrated data
% mesh is the calibrated mesh with initial guesses
% data link file will be reference for which data to use.
% mua_big and mus_big are the resulting values at each wavelength


if ~exist('nographs','var')
    nographs = 0;
end

parallel = parallel_init();

if parallel == 1
    nographs = 1;
end

% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

% If not a workspace variable, load meshes
if ischar(mesh_homog)== 1
  mesh_homog = load_mesh(mesh_homog);
end

if ischar(mesh_anom)== 1
  mesh_anom = load_mesh(mesh_anom);
end

mesh = mesh_anom;

% load anomaly data
paa_anom   = load_data(anom_data);

if ~exist('wv_array','var')
    wv_array = paa_anom.wv;
end

if ~isfield(paa_anom,'paa') || ~isfield(paa_anom,'wv')
    errordlg('Data not found or not properly formatted','NIRFAST Error');
    error('Data not found or not properly formatted');
end
if exist('wv_array','var')
    paa_anom_tmp.paa = [];
    paa_anom_tmp.link = paa_anom.link(:,1:2);
    paa_anom_tmp.wv = [];
    excoef_tmp = [];
    for i = 1:length(paa_anom.wv)
        if any(wv_array==paa_anom.wv(i))
            paa_anom_tmp.paa = [paa_anom_tmp.paa paa_anom.paa(:,i*2-1:i*2)];
            paa_anom_tmp.link = [paa_anom_tmp.link paa_anom.link(:,i+2)];
            paa_anom_tmp.wv = [paa_anom_tmp.wv paa_anom.wv(i)];
            a = mesh.wv == paa_anom.wv(i);
            excoef_tmp = [excoef_tmp; mesh.excoef(a,:)];
        end
    end
end
paa_anom = paa_anom_tmp; clear paa_anom_tmp;
mesh_anom.wv = paa_anom.wv;
mesh.wv = paa_anom.wv;
mesh.link = paa_anom.link;
mesh.excoef = excoef_tmp; clear excoef_tmp;
paa_anom = paa_anom.paa;

% Look for phase wrapping
[j,k] = size(paa_anom);
for i=1:2:k
paa_anom(paa_anom(:,i+1)<0,i+1) = ...
    paa_anom(paa_anom(:,i+1)<0,i+1) + (360);
paa_anom(paa_anom(:,i+1)>(360),i+1) = ...
    paa_anom(paa_anom(:,i+1)>(360),i+1) - (360);
end

% load homogeneous data
paa_homog  = load_data(homog_data);
if exist('wv_array','var')
    mesh_homog.wv = wv_array;
    paa_homog_tmp.paa = [];
    paa_homog_tmp.link = paa_homog.link(:,1:2);
    paa_homog_tmp.wv = [];
    excoef_tmp = [];
    for i = 1:length(paa_homog.wv)
        if any(wv_array==paa_homog.wv(i))
            paa_homog_tmp.paa = [paa_homog_tmp.paa paa_homog.paa(:,i*2-1:i*2)];
            paa_homog_tmp.link = [paa_homog_tmp.link paa_homog.link(:,i+2)];
            paa_homog_tmp.wv = [paa_homog_tmp.wv paa_homog.wv(i)];
            a = mesh_homog.wv == paa_homog.wv(i);
            excoef_tmp = [excoef_tmp; mesh_homog.excoef(a,:)];
        end
    end
end
paa_homog = paa_homog_tmp; clear paa_homog_tmp;
mesh_homog.link = paa_homog.link;
mesh_homog.wv = paa_homog.wv;
paa_homog = paa_homog.paa;
mesh_homog.excoef = excoef_tmp; clear excoef_tmp;

% Look for phase wrapping
[j,k] = size(paa_homog);
for i=1:2:k
paa_homog(paa_homog(:,i+1)<0,i+1) = ...
    paa_homog(paa_homog(:,i+1)<0,i+1) + (360);
paa_homog(paa_homog(:,i+1)>(360),i+1) = ...
    paa_homog(paa_homog(:,i+1)>(360),i+1) - (360);
end

mua_big = []; 
mus_big = [];
[j,k] = size(paa_homog);

% loop through wavelengths
kk = length(mesh.wv);
paa_homog_fit = reshape(paa_homog,j,2,kk);
paa_anom_fit = reshape(paa_anom,j,2,kk);

% SERIAL

if parallel == 0

    for i=1:kk
        % Make a standard mesh for this wavelength
        wvmesh = mesh_homog;
        wvmesh.type = 'stnd';
        wvmesh.link = [mesh_homog.link(:,1:2) mesh_homog.link(:,i+2)];
        [wvmesh.mua, wvmesh.mus, wvmesh.kappa] = calc_mua_mus(mesh_homog,mesh.wv(i));
        linkh = logical(wvmesh.link(:,3));
        % Calculate global mua and mus plus offsets for phantom data
        if frequency == 0
            paa_homog_fitcw = paa_homog_fit;
            [mua_h,mus_h,lnI_h,data_h_fem] = fit_data_cw(wvmesh,...
                                                              paa_homog_fitcw(:,1,i),...
                                                              iteration,...
                                                              nographs);
            phase_h = zeros(size(lnI_h,1),1);
        else
            [mua_h,mus_h,lnI_h,phase_h,data_h_fem] = fit_data(wvmesh,...
                                                              paa_homog_fit(:,:,i),...
                                                              frequency,...
                                                              iteration,...
                                                              nographs);
        end
       % mtit(['Homog ' num2str(mesh.wv(i)) 'nm'],'FontSize',14);
       disp(['Done homog ' num2str(mesh.wv(i)) 'nm'])

        wvmesh = [];

        % Make a standard mesh for this wavelength
        wvmesh = mesh_anom;
        wvmesh.type = 'stnd';
        wvmesh.link = [mesh.link(:,1:2) mesh.link(:,i+2)];
        [wvmesh.mua, wvmesh.mus, wvmesh.kappa] = calc_mua_mus(mesh,mesh.wv(i));
        linka = logical(wvmesh.link(:,3));
        % Calculate global mua and mus plus offsets for patient data
        if frequency == 0
            paa_anom_fitcw = paa_anom_fit;
            [mua_a,mus_a,lnI_a,data_a_fem] = fit_data_cw(wvmesh,...
                                                              paa_anom_fitcw(:,1,i),...
                                                              iteration,...
                                                              nographs);
            phase_a = zeros(size(lnI_a,1),1);
        else
            [mua_a,mus_a,lnI_a,phase_a,data_a_fem] = fit_data(wvmesh,...
                                                              paa_anom_fit(:,:,i),...
                                                              frequency,...
                                                              iteration,...
                                                              nographs);
        end
        % mtit(['Anom ' num2str(mesh.wv(i)) 'nm'],'FontSize',14); 
        disp(['Done anom ' num2str(mesh.wv(i)) 'nm'])

        wvmesh = [];



        % calculate offsets between modeled homogeneous and measured
        % homogeneous and using these calibrate data
        data_h_fem(:,1) = log(data_h_fem(:,1));
        paa_anomtmp = [log(paa_anom_fit(:,1,i)) paa_anom_fit(:,2,i)];
        paa_homogtmp = [log(paa_homog_fit(:,1,i)) paa_homog_fit(:,2,i)];

        paa_cal = paa_anomtmp - ((paa_homogtmp - data_h_fem));
        paa_cal(:,1) = paa_cal(:,1) - (lnI_a-lnI_h);
        paa_cal(:,2) = paa_cal(:,2) - (phase_a-phase_h);
        paa_cal(:,1) = exp(paa_cal(:,1));

        % calibrated data out into larger complete array
        data(:,:,i) = paa_cal;
        link(:,i) = and(mesh.link(:,i+2) , mesh_homog.link(:,i+2));

        % create an array for optical properties for all wavelengths
        mua_big = [mua_big; mua_a];
        mus_big = [mus_big; mus_a];

    end

else

% PARALLEL

    parfor i=1:kk
        % Make a standard mesh for this wavelength
        wvmesh = mesh_homog;
        wvmesh.type = 'stnd';
        wvmesh.link = [mesh_homog.link(:,1:2) mesh_homog.link(:,i+2)];
        [wvmesh.mua, wvmesh.mus, wvmesh.kappa] = calc_mua_mus(mesh_homog,mesh.wv(i));
        linkh = logical(wvmesh.link(:,3));
        % Calculate global mua and mus plus offsets for phantom data
        if frequency == 0
            paa_homog_fitcw = paa_homog_fit;
            [mua_h,mus_h,lnI_h,data_h_fem] = fit_data_cw(wvmesh,...
                                                              paa_homog_fitcw(:,1,i),...
                                                              iteration,...
                                                              nographs);
            phase_h = zeros(size(lnI_h,1),1);
        else
            [mua_h,mus_h,lnI_h,phase_h,data_h_fem] = fit_data(wvmesh,...
                                                              paa_homog_fit(:,:,i),...
                                                              frequency,...
                                                              iteration,...
                                                              nographs);
        end
       % mtit(['Homog ' num2str(mesh.wv(i)) 'nm'],'FontSize',14);
       disp(['Done homog ' num2str(mesh.wv(i)) 'nm'])

        wvmesh = [];

        % Make a standard mesh for this wavelength
        wvmesh = mesh_anom;
        wvmesh.type = 'stnd';
        wvmesh.link = [mesh.link(:,1:2) mesh.link(:,i+2)];
        [wvmesh.mua, wvmesh.mus, wvmesh.kappa] = calc_mua_mus(mesh,mesh.wv(i));
        linka = logical(wvmesh.link(:,3));
        % Calculate global mua and mus plus offsets for patient data
        if frequency == 0
            paa_anom_fitcw = paa_anom_fit;
            [mua_a,mus_a,lnI_a,data_a_fem] = fit_data_cw(wvmesh,...
                                                              paa_anom_fitcw(:,1,i),...
                                                              iteration,...
                                                              nographs);
            phase_a = zeros(size(lnI_a,1),1);
        else
            [mua_a,mus_a,lnI_a,phase_a,data_a_fem] = fit_data(wvmesh,...
                                                              paa_anom_fit(:,:,i),...
                                                              frequency,...
                                                              iteration,...
                                                              nographs);
        end
        % mtit(['Anom ' num2str(mesh.wv(i)) 'nm'],'FontSize',14); 
        disp(['Done anom ' num2str(mesh.wv(i)) 'nm'])

        wvmesh = [];



        % calculate offsets between modeled homogeneous and measured
        % homogeneous and using these calibrate data
        data_h_fem(:,1) = log(data_h_fem(:,1));
        paa_anomtmp = [log(paa_anom_fit(:,1,i)) paa_anom_fit(:,2,i)];
        paa_homogtmp = [log(paa_homog_fit(:,1,i)) paa_homog_fit(:,2,i)];

        paa_cal = paa_anomtmp - ((paa_homogtmp - data_h_fem));
        paa_cal(:,1) = paa_cal(:,1) - (lnI_a-lnI_h);
        paa_cal(:,2) = paa_cal(:,2) - (phase_a-phase_h);
        paa_cal(:,1) = exp(paa_cal(:,1));

        % calibrated data out into larger complete array
        data(:,:,i) = paa_cal;
        link(:,i) = and(mesh.link(:,i+2) , mesh_homog.link(:,i+2));

        % create an array for optical properties for all wavelengths
        mua_big = [mua_big; mua_a];
        mus_big = [mus_big; mus_a];

    end

end

data.paa = reshape(data,j,k);
data.link = [mesh.link(:,1:2) link];

% set wavelengths array
data.wv = mesh.wv;
wv_array = data.wv';
mesh.link = data.link;

[nc,junk] = size(mesh.chromscattlist);
nwn = length(wv_array);

% get extinction coeffs for chosen wavelengths
[junk, junk1, junk2, E] = calc_mua_mus(mesh,wv_array);
clear junk junk1 junk2;

% using a non-negative least squares fit, get concentrations.
C = lsqnonneg(E,mua_big);
mesh.conc = repmat(C',length(mesh.nodes),1);

% now fit for particle size (a) and density (b)
% note: mus = a* lambda^-b
%       log(mus) = log(a) - b*(log(lambda))
% so fit for gradient (-b) and offset (log(a))
xdata = wv_array./1000;
ydata = mus_big;
p = polyfit(log(xdata),log(ydata),1);

mesh.sp(:,1) = -p(1);
mesh.sa(:,1) = exp(p(2));

% generate initial guess if optimization toolbox exists
if exist('constrain_lsqfit') % it's fit...
     display('Using Optimization Toolbox');
    [A,b,Aeq,beq,lb,ub]=constrain_lsqfit(nwn,nc);

    [C,resnorm,residual,exitflag,output,lambda] = lsqlin(E,mua_big,A,b,Aeq,beq,lb,ub);


    if(C(end)==0)
        % for water recon (no fat)...
        ub(end) = 0.5;lb(end) = 0.5;
        [C,resnorm,residual,exitflag,output,lambda] = lsqlin(E,mua_big,A,b,Aeq,beq,lb,ub);
    end
    
    mesh.conc = repmat(C',length(mesh.nodes),1);

    % for scattering parameters
    x0 = [1.0;1.0];
    xdata = wv_array./1000;

    [x,resnorm] = lsqcurvefit(@power_fit,x0,xdata,mus_big,[.2 .2],[3 3])

    mesh.sa(:) = x(1);
    mesh.sp(:) = x(2);

end
