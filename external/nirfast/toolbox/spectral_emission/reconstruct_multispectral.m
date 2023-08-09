function [mesh,pj_error_fl] = reconstruct_multispectral(mesh,recon_basis,frequency,data_fn,...
            iteration,lambda,output_fn,filter_n,drug,wv_array_emiss,wv_excite,eta_tot_true,tau)
        
% [mesh,pj_error_fl] = reconstruct_multispectral(mesh,recon_basis,frequency,data_fn,...
%     iteration,lambda,output_fn,filter_n,drug,wv_array_emiss,wv_excite,eta_tot_true,tau)
%
% spectral emission reconstruction, find fluorescence info based on
% boundary data
%
% mesh is the input mesh (variable or filename)
% recon_basis is the reconstruction basis (pixel basis or mesh filename)
% frequency is the modulation frequency (MHz)
% data_fn is the boundary data (variable or filename)
% iteration is the max number of iterations
% lambda is the initial regularization value
% output_fn is the root output filename
% filter_n is the number of mean filters
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
if ischar(mesh)== 1
  mesh = load_mesh(mesh);
end
if ~strcmp(mesh.type,'spec')
    errordlg('Mesh type must be spectral (''spec'')','NIRFAST Error');
    error('Mesh type must be spectral (''spec'')');
end

%% Initiate log file
etamuaf_sol=[output_fn '_etamuaf.sol'];
fid_log = fopen([output_fn '.log'],'w');
fprintf(fid_log,'Forward Mesh   = %s\n',mesh.name);
if ischar(recon_basis)
    fprintf(fid_log,'Basis          = %s\n',recon_basis);
end
fprintf(fid_log,'Frequency      = %f MHz\n',frequency);
if ischar(data_fn) ~= 0
    fprintf(fid_log,'Data File      = %s\n',data_fn);
end
if isstruct(lambda)
    fprintf(fid_log,'Initial Regularization  = %d\n',lambda.value);
else
    fprintf(fid_log,'Initial Regularization  = %d\n',lambda);
end
fprintf(fid_log,'Filtering        = %d\n',filter_n);
fprintf(fid_log,'Output Files   = %s',etamuaf_sol);
fprintf(fid_log,'Initial Guess muaf = %d\n',mesh.muaf(1));
fprintf(fid_log,'\n');

%% load recon_mesh
disp('Loading recon basis')
if ischar(recon_basis)
  recon_mesh = load_mesh(recon_basis);
  [mesh.fine2coarse,...
   recon_mesh.coarse2fine] = second_mesh_basis(mesh,recon_mesh);
elseif isstruct(recon_basis)
  recon_mesh = recon_basis;
  [mesh.fine2coarse,...
   recon_mesh.coarse2fine] = second_mesh_basis(mesh,recon_mesh);
else
  [mesh.fine2coarse,recon_mesh] = pixel_basis(recon_basis,mesh);
end

%% load data
anom = load_data(data_fn,wv_array_emiss);
if isempty(anom) || ~isfield(anom,'paa')
    errordlg('Data not found or not properly formatted','NIRFAST Error');
    error('Data not found or not properly formatted');
end

% match link 0s with data NaNs
datanum = 0;
[ns,junk]=size(mesh.source.coord);
for i = 1 : ns
    for j = 1 : length(mesh.link(i,:))
        datanum = datanum+1;
        if sum(isnan(anom.paa(datanum,:)))>0
            mesh.link(i,j) = 0;
        end
        if mesh.link(i,j) == 0
            anom.paa(datanum,:) = NaN;
        end
    end
end

% we need log amplitude
[nr,nc]=size(anom.paa);
k = 1;
for i = 1 : 2 : nc
    anom_a(:,k) = log(anom.paa(:,i));
    k = k + 1;
end

% find NaN in data
for i = 1 : length(wv_array_emiss)
    tmp = find(isnan(anom_a(:,i))==1);
    eval(['ind.l' num2str(wv_array_emiss(i)) ' = tmp;']);
end
clear tmp

% Create a vector of data
anom = [];
for i = 1 : length(wv_array_emiss)
    eval(['tmp = ind.l' num2str(wv_array_emiss(i)) ';']);
    ind_tmp = setdiff(1:size(anom_a(:,i),1),tmp);
    tmp = [anom_a(ind_tmp,i)];
    [nr,nc]=size(tmp);
    anom = [anom; reshape(tmp',nc*nr,1)];
end
mesh.ind = ind;
clear anom_* tmp ind*

%% load emission spectrum
emiss_spec_temp = load(drug);
emiss_spec.wv = emiss_spec_temp(:,1);
emiss_spec.etaspec = emiss_spec_temp(:,2);
clear emiss_spec_temp;
 
%% recon parameters

% calculate initial guess of muaf
ex = mesh.excoef(find(mesh.wv == wv_excite),end);
mesh.muaf = mesh.conc(:,end)*ex;

% Currently, all fluorescence data is in CW 
omega = 2*pi*frequency*1e6;

% calculate excitation parameters
[mesh.muax, mesh.musx, mesh.kappax, junk] = calc_mua_mus(mesh,wv_excite);

% set tau
mesh.tau = zeros(size(mesh.muax,1),1);
mesh.tau(:) = tau;

% initialize projection error
pj_error_fl=[];

%% excitation data
mesh.mua = mesh.muax;
mesh.mus = mesh.musx;
mesh.kappa = mesh.kappax;
datax = femdata_stnd(mesh,100);

%% check for input regularization
if isstruct(lambda) && ~(strcmp(lambda.type,'JJt') || strcmp(lambda.type,'JtJ'))
    lambda.type = 'Automatic';
end
if ~isstruct(lambda)
    lambda.value = lambda;
    lambda.type = 'Automatic';
end
% determine regularization type
if strcmp(lambda.type, 'Automatic')
    if size(anom,1)<2*size(recon_mesh.nodes,1)
        lambda.type = 'JJt';
    else
        lambda.type = 'JtJ';
    end
end

%% iterate
for it=1:iteration

    %% calculate jacobian
    Jbig = [];
    Ref = [];
    for i = 1:numel(wv_array_emiss)

        disp(sprintf('Calculate Jacobian for emission field at %g nm',wv_array_emiss(i))); 

        % determine mua mus values from concetrations and scatt params
        [mesh.muam, mesh.musm, mesh.kappam, junk] = calc_mua_mus(mesh,wv_array_emiss(i));
        mesh.eta = zeros(size(mesh.muam,1),1);
        mesh.eta(:) = emiss_spec.etaspec(find(emiss_spec.wv == wv_array_emiss(i)))*eta_tot_true;

        % set fluorescence variables
        mesh.gamma = (mesh.eta.*mesh.muaf)./(1+(omega.*mesh.tau).^2);

        % build jacobian
        [J,datafl] = jacobian_fl(mesh,frequency,datax);
        [nrow1,junk] = size(J.completem);
        
        % extract ref data
        ref(:,1) = log(datafl.amplitudem);
        Ref = [Ref; ref];
        
        % Interpolate Jacobian onto recon mesh
        [jm,recon_mesh] = interpolatef2r_fl(mesh,recon_mesh,J.completem);
        jm = jm(1:2:end-1, 1:end/2); % take only intensity portion

        % Normalize Jacobian wrt fl source gamma
        jm = jm*diag([recon_mesh.gamma]);
        
        Jbig = [Jbig; jm];
        clear ref J
    end
      
    % Screen and Log Info
    pj_error_fl = [pj_error_fl sum((anom-Ref).^2)];
    
    disp('---------------------------------');
    disp(['Iteration_fl Number          = ' num2str(it)]);
    disp(['Projection_fl error          = ' num2str(pj_error_fl(end))]);
    
    fprintf(fid_log,'---------------------------------\n');
    fprintf(fid_log,'Iteration_fl Number          = %d\n',it);
    fprintf(fid_log,'Projection_fl error          = %f\n',pj_error_fl(end));
    
    if it ~= 1
        p = (pj_error_fl(end-1)-pj_error_fl(end))*100/pj_error_fl(end-1);
        disp(['Projection_fl error change   = ' num2str(p) '%']);
        fprintf(fid_log,'Projection_fl error change   = %f %%\n',p);
        if (p) <= 2
            disp('---------------------------------');
            disp('STOPPING CRITERIA FOR FLUORESCENCE COMPONENT REACHED');
            fprintf(fid_log,'---------------------------------\n');
            fprintf(fid_log,'STOPPING CRITERIA FOR FLUORESCENCE COMPONENT REACHED\n');
            % set output
            data_recon.elements = mesh.elements;
            data_recon.etamuaf = mesh.etamuaf;  
            break
        end
    end
    clear data_recon

    if strcmp(lambda.type, 'JJt')
        %% build Hessian
        [nrow,ncol]=size(Jbig);
        Hess = zeros(ncol);
        Hess = Jbig*Jbig';

        % add regularization
        reg = lambda.value.*(max(diag(Hess)));
        disp(['Regularization Fluor           = ' num2str(reg)]);
        fprintf(fid_log,'Regularization Fluor            = %f\n',reg);
        Hess = Hess+(eye(nrow).*reg);

        %% Calculate update
        data_diff = (anom-Ref);
        u = Jbig'*(Hess\data_diff);
        u = u.*[recon_mesh.gamma];
    else
        %% build Hessian
        [nrow,ncol]=size(Jbig);
        Hess = zeros(ncol);
        Hess = Jbig'*Jbig;

        % add regularization
        reg = lambda.value.*(max(diag(Hess)));
        disp(['Regularization Fluor           = ' num2str(reg)]);
        fprintf(fid_log,'Regularization Fluor            = %f\n',reg);
        Hess = Hess+(eye(ncol).*reg);

        %% Calculate update
        data_diff = (anom-Ref);
        u = (Hess\Jbig'*data_diff);
        u = u.*[recon_mesh.gamma];
    end

    % value update:  
    recon_mesh.gamma = recon_mesh.gamma+u;
    recon_mesh.etamuaf = recon_mesh.gamma.*(1+(omega.*recon_mesh.tau).^2);
    recon_mesh.etamuaf(find(recon_mesh.etamuaf < 0)) = 0;
    recon_mesh.muaf = recon_mesh.etamuaf./recon_mesh.eta;
    clear u Hess Hess_norm tmp data_diff G
    
    % interpolate onto fine mesh  
    [mesh,recon_mesh] = interpolatep2f_fl(mesh,recon_mesh);
    
    % filter
    if filter_n ~= 0
        disp('Filtering');
        fwd_mesh = mean_filter(mesh,filter_n);
    end
    
    % Write solution to file
    if it == 1
        fid = fopen(etamuaf_sol,'w');
    else
        fid = fopen(etamuaf_sol,'a');
    end
    fprintf(fid,'solution %d ',it);
    fprintf(fid,'-size=%g ',length(mesh.nodes));
    fprintf(fid,'-components=1 ');
    fprintf(fid,'-type=nodal\n');
    fprintf(fid,'%g ',mesh.etamuaf);
    fprintf(fid,'\n');
    fclose(fid);
end

fin_it = it-1;
mesh.type='fluor';



%% Sub functions
function [val_int,recon_mesh] = interpolatef2r_fl(fwd_mesh,recon_mesh,val)

% This function interpolates fwd_mesh into recon_mesh
% For the Jacobian it is an integration!
NNC = size(recon_mesh.nodes,1);
NNF = size(fwd_mesh.nodes,1);
NROW = size(val,1);
val_int = zeros(NROW,NNC*2);

for i = 1 : NNF
    if recon_mesh.coarse2fine(i,1) ~= 0
        val_int(:,recon_mesh.elements(recon_mesh.coarse2fine(i,1),:)) = ...
            val_int(:,recon_mesh.elements(recon_mesh.coarse2fine(i,1),:)) + ...
            val(:,i)*recon_mesh.coarse2fine(i,2:end);
        val_int(:,recon_mesh.elements(recon_mesh.coarse2fine(i,1),:)+NNC) = ...
            val_int(:,recon_mesh.elements(recon_mesh.coarse2fine(i,1),:)+NNC) + ...
            val(:,i+NNF)*recon_mesh.coarse2fine(i,2:end);
    elseif recon_mesh.coarse2fine(i,1) == 0
        dist = distance(fwd_mesh.nodes,fwd_mesh.bndvtx,recon_mesh.nodes(i,:));
        mindist = find(dist==min(dist));
        mindist = mindist(1);
        val_int(:,i) = val(:,mindist);
        val_int(:,i+NNC) = val(:,mindist+NNF);
    end
end

for i = 1 : NNC
    if fwd_mesh.fine2coarse(i,1) ~= 0
        recon_mesh.region(i,1) = ...
            median(fwd_mesh.region(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:)));
        recon_mesh.eta(i,1) = (fwd_mesh.fine2coarse(i,2:end) * ...
            fwd_mesh.eta(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:)));
        recon_mesh.muaf(i,1) = (fwd_mesh.fine2coarse(i,2:end) * ...
            fwd_mesh.muaf(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:)));
        recon_mesh.gamma(i,1) = (fwd_mesh.fine2coarse(i,2:end) * ...
            fwd_mesh.gamma(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:)));
        recon_mesh.tau(i,1) = (fwd_mesh.fine2coarse(i,2:end) * ...
            fwd_mesh.tau(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:)));
        
    elseif fwd_mesh.fine2coarse(i,1) == 0
        dist = distance(fwd_mesh.nodes,...
            fwd_mesh.bndvtx,...
            [recon_mesh.nodes(i,1:2) 0]);
        mindist = find(dist==min(dist));
        mindist = mindist(1);
        recon_mesh.region(i,1) = fwd_mesh.region(mindist);
        recon_mesh.eta(i,1) = fwd_mesh.eta(mindist);
        recon_mesh.muaf(i,1) = fwd_mesh.muaf(mindist);
        recon_mesh.gamma(i,1) = fwd_mesh.gamma(mindist);
        recon_mesh.tau(i,1) = fwd_mesh.tau(mindist);
        
    end
end


function [fwd_mesh,recon_mesh] = interpolatep2f_fl(fwd_mesh,recon_mesh)

for i = 1 : length(fwd_mesh.nodes)
  fwd_mesh.gamma(i,1) = ...
      (recon_mesh.coarse2fine(i,2:end) * ...
       recon_mesh.gamma(recon_mesh.elements(recon_mesh.coarse2fine(i,1),:)));
  fwd_mesh.muaf(i,1) = ...
      (recon_mesh.coarse2fine(i,2:end) * ...
       recon_mesh.muaf(recon_mesh.elements(recon_mesh.coarse2fine(i,1),:)));
   fwd_mesh.eta(i,1) = ...
      (recon_mesh.coarse2fine(i,2:end) * ...
       recon_mesh.eta(recon_mesh.elements(recon_mesh.coarse2fine(i,1),:)));
  fwd_mesh.etamuaf(i,1) = ...
      (recon_mesh.coarse2fine(i,2:end) * ...
       recon_mesh.etamuaf(recon_mesh.elements(recon_mesh.coarse2fine(i,1),:)));
  fwd_mesh.tau(i,1) = ...
      (recon_mesh.coarse2fine(i,2:end) * ...
       recon_mesh.tau(recon_mesh.elements(recon_mesh.coarse2fine(i,1),:)));
end