function [fwd_mesh,pj_error] = reconstruct_spectral_spatial(fwd_mesh,...
    recon_basis,...
    frequency,...
    data_fn,...
    iteration,...
    lambda,...
    output_fn,...
    filter_n,...
    wv_array)

% [fwd_mesh,pj_error] = reconstruct_spectral(fwd_mesh,...
%                                            recon_basis,...
%                                            frequency,...
%                                            data_fn,...
%                                            iteration,...
%                                            lambda,...
%                                            output_fn,...
%                                            filter_n,...
%                                            wv_array)
%
% Spectrally constrained or "direct spectral" reconstruction program.
%
% fwd_mesh is the input mesh (variable or filename)
% recon_basis is the reconstruction basis (pixel or mesh)
% frequency is the modulation frequency (MHz)
% data_fn is the data filename
% iteration is the max number of iterations
% lambda is the initial regularization value
% output_fn is the root output filename
% filter_n is the number of mean filters
% wv_array is optional wavelength array


tic

% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

%*******************************************************
% If not a workspace variable, load mesh
if ischar(fwd_mesh)== 1
    fwd_mesh = load_mesh(fwd_mesh);
end
if ~strcmp(fwd_mesh.type,'spec')
    errordlg('Mesh type is incorrect','NIRFAST Error');
    error('Mesh type is incorrect');
end
if exist('wv_array') == 0
    wv_array = fwd_mesh.wv;
end

% check to ensure wv_array wavelengths match the wavelength list fwd_mesh
for i = 1:length(wv_array)
    tmp = find(fwd_mesh.wv == wv_array(i));
    if isempty(tmp)
        flag(i) = 0;
    else
        flag(i) = tmp(1);
    end
end
tmp = find(flag==0);
if isempty(tmp) ~= 1
    for i = 1 : length(tmp)
        disp(['ERROR: wv_array contains ' num2str(wv_array(tmp(i))) ...
            'nm which is not present in ' fwd_mesh.name,'.excoef']);
    end
    return
end
clear tmp flag i
nwv = length(wv_array);

%*******************************************************
% read data - This is the calibrated experimental data or simulated data
disp('Loading data and wavelength information')
data = load_data(data_fn,wv_array);
% if specified wavelength not available, terminate.
if isempty(data) || ~isfield(data,'paa')
    errordlg('Data not found or not properly formatted','NIRFAST Error');
    error('Data not found or not properly formatted');
end
[n,m] = size(data.link);  [nd,md] = size(data.paa);
data_link = data.link;

% we need log amplitude and phase in radians
anom_a = [];
anom_p = [];
k=1;
for i = 1:2:md
    data.paa(:,i) = log(data.paa(:,i));
    foo = data.paa(:,i+1)/180.0*pi;
    foo(foo<0) = foo(foo<0) + (2*pi);
    foo(foo>(2*pi)) = foo(foo>(2*pi)) - (2*pi);
    data.paa(:,i+1) = foo; clear foo
    linki = logical(data_link(:,k+2));
    anom_a = [anom_a; data.paa(linki,i)];
    anom_p = [anom_p; data.paa(linki,i+1)];
    k = k+1;
end
anom = zeros(length(anom_a)*2,1);
anom(1:2:end) = anom_a;
anom(2:2:end) = anom_p;
clear data anom_a anom_p

% extinction coeff for chosen wavelengths
[junk1,junk2,junk3,E] = calc_mua_mus(fwd_mesh,wv_array);
clear junk*

fwd_mesh.link = data_link;
%*******************************************************
% initialize projection error
pj_error = [];

%*******************************************************
% load recon_mesh
disp('Loading recon basis')
if ischar(recon_basis)
    recon_mesh = load_mesh(recon_basis);
    [fwd_mesh.fine2coarse,...
        recon_mesh.coarse2fine] = second_mesh_basis(fwd_mesh,recon_mesh);
elseif isstruct(recon_basis)
    recon_mesh = recon_basis;
    [fwd_mesh.fine2coarse,...
        recon_mesh.coarse2fine] = second_mesh_basis(fwd_mesh,recon_mesh);
else
    [fwd_mesh.fine2coarse,recon_mesh] = pixel_basis(recon_basis,fwd_mesh);
end

% set parameters for second mesh. Jacobian calculation has been modified so
% that it only calculates Jacobian for second mesh, which is more memory
% efficient as well as computationally faster
recon_mesh.type = fwd_mesh.type;
recon_mesh.link = fwd_mesh.link;
recon_mesh.source = fwd_mesh.source;
recon_mesh.meas = fwd_mesh.meas;
recon_mesh.dimension = fwd_mesh.dimension;
recon_mesh.wv = fwd_mesh.wv;
recon_mesh.excoef = fwd_mesh.excoef;
if recon_mesh.dimension == 2
    recon_mesh.element_area = ele_area_c(recon_mesh.nodes(:,1:2),...
        recon_mesh.elements);
else
    recon_mesh.element_area = ele_area_c(recon_mesh.nodes,...
        recon_mesh.elements);
end
pixel.support = mesh_support(recon_mesh.nodes,...
    recon_mesh.elements,...
    recon_mesh.element_area);

%************************************************
% Initiate log file
fid_log = fopen([output_fn,'.log'],'w');
fprintf(fid_log,'Forward Mesh       = %s\n',fwd_mesh.name);
if ischar(recon_basis)
    fprintf(fid_log,'Basis              = %s\n',recon_basis);
end
fprintf(fid_log,'Frequency          = %f MHz\n',frequency);
if ischar(data_fn) ~= 0
    fprintf(fid_log,'Data file          = %s\n',data_fn);
end
if isstruct(lambda)
    fprintf(fid_log,'Initial Regularization  = %d\n',lambda.value);
else
    fprintf(fid_log,'Initial Regularization  = %d\n',lambda);
end
fprintf(fid_log,'Filter             = %d\n',filter_n);
fprintf(fid_log,'Wavelengths Used   = ');
for i = 1 : nwv
    fprintf(fid_log,'%s ',num2str(fwd_mesh.wv(i)));
end
fprintf(fid_log,'\n\n');
fprintf(fid_log,'********************************\n');
fprintf(fid_log,'Output files:\n');

% Get string names for Chromophore files
all_sol = char(fwd_mesh.chromscattlist);
[n_allsol,junk] = size(all_sol);
str_label = [];
for i = 1:n_allsol
    tmp_string = [strcat(output_fn,'_',all_sol(i,:),'.sol')];
    str_label = [str_label tmp_string];
    fprintf(fid_log,'%s\n',tmp_string);
    clear tmp_string
end

% log initial guesses
for i = 1:n_allsol
    if strcmp(all_sol(i,1:3),'S-A')
        fprintf(fid_log,['Initial Guess ' all_sol(i,:) ' = %d\n'],...
            fwd_mesh.sa(1));
    elseif strcmp(all_sol(i,1:3),'S-P')
        fprintf(fid_log,['Initial Guess ' all_sol(i,:) ' = %d\n'],...
            fwd_mesh.sp(1));
    else
        fprintf(fid_log,['Initial Guess ' all_sol(i,:) ' = %d\n'],...
            fwd_mesh.conc(1,i));
    end
end


% start non-linear itertaion image reconstruction part
for it = 1:iteration
    
    % Interpolate fwd_mesh values to recon_mesh
    [recon_mesh.conc,recon_mesh.sa,recon_mesh.sp] = ...
        interp2coarse(fwd_mesh,...
        recon_mesh,...
        fwd_mesh.conc,...
        fwd_mesh.sa,...
        fwd_mesh.sp);
    
    % Compute Jacobian and ref data from each wavelength data
    disp('---------------------------------');
    disp('Building Jacobian using jacobian_spectral')
    [J,data,fwd_mesh] = jacobian_spectral(fwd_mesh,frequency,wv_array,recon_mesh);
    
    % we need log amplitude and phase in radians    
    ref_a = [];
    ref_p = [];
    k=1;
    for i = 1:2:md
        data.paa(:,i) = log(data.paa(:,i));
        foo = data.paa(:,i+1)/180.0*pi;
        foo(foo<0) = foo(foo<0) + (2*pi);
        foo(foo>(2*pi)) = foo(foo>(2*pi)) - (2*pi);
        data.paa(:,i+1) = foo; clear foo
        linki = logical(data_link(:,k+2));
        ref_a = [ref_a; data.paa(linki,i)];
        ref_p = [ref_p; data.paa(linki,i+1)];
        k = k+1;
    end
    ref = zeros(length(ref_a)*2,1);
    ref(1:2:end) = ref_a;
    ref(2:2:end) = ref_p;
    clear data ref_a ref_p

    
    % Calculate data difference:
    data_diff = (anom - ref);
    
    %*********************************************************************
    % Update pj_error and check stopping criteria
    pj_error = [pj_error sum(abs(data_diff.^2))];
    
    disp(['Iteration Number          = ' num2str(it)]);
    disp(['Projection error          = ' num2str(pj_error(end))]);
    
    fprintf(fid_log,'---------------------------------\n');
    fprintf(fid_log,'Iteration Number          = %d\n',it);
    fprintf(fid_log,'Projection error          = %f\n',pj_error(end));
    
    if it ~= 1
        p = (pj_error(end-1)-pj_error(end))*100/pj_error(end-1);
        disp(['Projection error change   = ' num2str(p) '%']);
        fprintf(fid_log,'Projection error change   = %f %%\n',p);
        if (p) <= 2
            disp('---------------------------------');
            disp('STOPPING CRITERIA REACHED');
            fprintf(fid_log,'---------------------------------\n');
            fprintf(fid_log,'STOPPING CRITERIA REACHED\n');
            break
        end
    end

    
    %*********************************************************************
    % Normalize Jacobian w.r.t different SI units
  conc_temp = reshape(recon_mesh.conc, numel(recon_mesh.conc), 1);
   N = [conc_temp; recon_mesh.sa; recon_mesh.sp];
  [nr,nc]=size(J);
  for i = 1 : nr
    J(i,:) = J(i,:).*N';
  end
    
    % build hessian
    [nrow,ncol]=size(J);
        
    
        Hess = zeros(ncol);
        disp('Calculating Hessian');
        Hess = (J'*J);
        
 % Calculate L-matrix
  if it == 1
    disp('calculating L matrix');
    tic
    L = create_L_fast_final(length(recon_mesh.nodes),...
			    length(unique(recon_mesh.region)),...
			    recon_mesh.region);
    LtL = L'*L;
    clear L L_paa
    toc
  end
  
  % Add regularization 
  nn = length(recon_mesh.nodes);
  nc = length(Hess)/nn;
  reg = [];
  for i = 1 : nc
    reg(i) = lambda.*max(diag(Hess((i-1)*nn+1:i*nn,(i-1)*nn+1:i*nn)));
  end
 
  disp(['Regularization        = ' num2str(reg)]);
  fprintf(fid_log,'Regularization        = %s\n',num2str(reg));

  for i = 1 : nc
    reg_LtL = LtL.*reg(i);
    Hess((i-1)*nn+1:i*nn,(i-1)*nn+1:i*nn) = ...
	Hess((i-1)*nn+1:i*nn,(i-1)*nn+1:i*nn) + reg_LtL;
  end
  
  % Calculate update
  %foo = minres(Hess,J'*data_diff,1e-5,2000);
  foo = Hess\(J'*(data_diff));
  clear J reg Hess;
    
    %******************************************************
    % Inversion complete, normalize update
     foo = foo.*N;
     
    % Update values
    [nn,nc] = size(recon_mesh.conc);
    foo_conc = foo(1:nn*nc);
    foo_sa = foo(nn*nc+1:nn*(nc+1));
    foo_sp = foo(nn*(nc+1)+1:end);
    
    conc_temp = conc_temp + foo_conc;
    recon_mesh.conc = reshape(conc_temp,nn,nc);
    recon_mesh.sa = recon_mesh.sa + foo_sa;
    recon_mesh.sp = recon_mesh.sp + foo_sp;
    clear foo foo_conc foo_sa foo_sp conc_temp N
    
    %% constraining water to be less than 100%, sa and sp <3.0
    [recon_mesh.conc, recon_mesh.sa, recon_mesh.sp] = ...
        constrain_val(recon_mesh, recon_mesh.conc, ...
        recon_mesh.sa, recon_mesh.sp, ...
        fwd_mesh.chromscattlist);
    
    % Interpolate updated values back to fwd_mesh
    [fwd_mesh.conc, fwd_mesh.sa, fwd_mesh.sp] = ...
        interp2fine(fwd_mesh,...
        recon_mesh,...
        recon_mesh.conc,...
        recon_mesh.sa,...
        recon_mesh.sp);
    
    %%filtering
    if (filter_n ~= 0)
        disp('Filtering');
        [fwd_mesh] = mean_filter_chromscatt(fwd_mesh, filter_n);
    end
    
    
    %**********************************************************
    % Compute absorption and scatter coefficients
    % mua From Beer's Law:
    fwd_mesh.mua = (E*fwd_mesh.conc')';
    
    % mus' From Mie Theory:
    for k =1 : nwv
        fwd_mesh.mus(:,k) = ((fwd_mesh.sa).*(wv_array(k)/1000).^(-fwd_mesh.sp));
    end
    
    % Diffusion coeff:
    fwd_mesh.kappa = (1./(3*(fwd_mesh.mus + fwd_mesh.mua)));
    
    %****************************************************************
    % Calculate Clinical Parameters
    chrom_scatt = fwd_mesh.conc;
    chrom_scatt(:,end+1) = fwd_mesh.sa;
    chrom_scatt(:,end+1) = fwd_mesh.sp;
    
    %Hbt and StO2
    %chrom_scatt(:,1) = (fwd_mesh.conc(:,1) + fwd_mesh.conc(:,2)).*1000;
    %chrom_scatt(:,2) = ((fwd_mesh.conc(:,1).*1000)./chrom_scatt(:,1))*100;
    
    %****************************************************************
    % Save solutions to file
    
    for j = 1:n_allsol
        tmp_string = [strcat(output_fn,'_',all_sol(j,:),'.sol')];
        if (it==1)
            fid = fopen(tmp_string,'w');
        else
            fid = fopen(tmp_string,'a');
        end
        fprintf(fid,'solution %d ',it);
        fprintf(fid,'-size=%g ',length(fwd_mesh.nodes));
        fprintf(fid,'-components=1 ');
        fprintf(fid,'-type=nodal\n');
        fprintf(fid,'%f ',chrom_scatt(:,j));
        fprintf(fid,'\n');
        fclose(fid);
        clear tmp_string
    end
    clear chrom_scatt
    pause(1);
    plotmesh(fwd_mesh);
end
% Close log file
time = toc;
fprintf(fid_log,'Computation Time = %f\n',time);
fprintf(fid_log,'Final Solution:\n');
fprintf(fid_log,'%g',it-1);
fclose(fid_log);

% END RECONSTRUCTION
%**************************************************************************




%**************************************************************************
% Selected Sub-functions

function [conc_coarse,sa_coarse,sp_coarse] = interp2coarse(mesh1,mesh2,conc_f,sa_f,sp_f)
[junk,m] = size(conc_f);
for  i = 1:length(mesh2.nodes)
    if mesh1.fine2coarse(i,1) ~= 0
        for j = 1:m
            conc_coarse(i,j) = (mesh1.fine2coarse(i,2:end) * ...
                conc_f(mesh1.elements(mesh1.fine2coarse(i,1),:),j));
        end
        sa_coarse(i,1) = (mesh1.fine2coarse(i,2:end) * ...
            sa_f(mesh1.elements(mesh1.fine2coarse(i,1),:)));
        sp_coarse(i,1) = (mesh1.fine2coarse(i,2:end) * ...
            sp_f(mesh1.elements(mesh1.fine2coarse(i,1),:)));
        
    elseif mesh1.fine2coarse(i,1) == 0
        dist = distance(mesh1.nodes,...
            mesh1.bndvtx,...
            [mesh2.nodes(i,1:2) 0]);
        mindist = find(dist==min(dist));
        mindist = mindist(1);
        for j = 1:m
            conc_coarse(i,j) = conc_f(mindist,j);
        end
        sa_coarse(i,1) = sa_f(mindist);
        sp_coarse(i,1) = sp_f(mindist);
    end
end


function [conc_f,sa_f,sp_f] = interp2fine(mesh1,mesh2,conc_coarse,sa_coarse,sp_coarse)
[junk,m] = size(conc_coarse);
for  i = 1:length(mesh1.nodes)
    for j = 1:m
        conc_f(i,j) = (mesh2.coarse2fine(i,2:end) * ...
            conc_coarse(mesh2.elements(mesh2.coarse2fine(i,1),:),j));
    end
    sa_f(i,1) = (mesh2.coarse2fine(i,2:end) * ...
        sa_coarse(mesh2.elements(mesh2.coarse2fine(i,1),:)));
    sp_f(i,1) = (mesh2.coarse2fine(i,2:end) * ...
        sp_coarse(mesh2.elements(mesh2.coarse2fine(i,1),:)));
end


function [conc,sa,sp] = constrain_val(mesh2,conc,sa,sp,list)
% Constrain water
list = char(list);
list = list(:,1:5);
[nr,nc]=size(list);
for i = 1 : nr-2
    if strcmp(list(i,:),'Water') == 1
        index = find(conc(:,i) > 1.0);
        conc(index,i) = 1.0;
        clear index;
        index = find(conc(:,i) < 0.0);
        conc(index,i) = 0.0001;
        clear index;
    else
        index = find(conc(:,i) < 0.0);
        conc(index,i) = 0.0001;
        clear index;
    end
end

%%constraining scatt ampl
index = find(sa > 3.0);
sa(index) = 3.0;
clear index;
index = find(sa < 0.0);
sa(index) = 0.0001;
clear index;

%%constraining scatt power
index = find(sp > 3.0);
sp(index) = 3.0;
clear index;
index = find(sp < 0.0);
sp(index) = 0.0001;
clear index;
