function [fwd_mesh,pj_error] = reconstruct_fl(fwd_mesh,...
    recon_basis,...
    frequency,...
    data_fn,...
    iteration,...
    lambda,...
    output_fn,...
    filter_n)

% [fwd_mesh,pj_error] = reconstruct_fl(fwd_mesh,...
%                                      recon_basis,...
%                                      frequency,...
%                                      data_fn,...
%                                      iteration,...
%                                      lambda,...
%                                      output_fn,...
%                                      filter_n)
%
% reconstruction program for fluorescence meshes
%
% fwd_mesh is the input mesh (variable or filename)
% recon_basis is the reconstruction basis (pixel basis or mesh filename)
% frequency is the modulation frequency (MHz)
% data_fn is the boundary data (variable or filename)
% iteration is the max number of iterations
% lambda is the initial regularization value
% output_fn is the root output filename
% filter_n is the number of mean filters



% always CW for fluor
frequency = 0;

%*******************************************************
% Read data
data = load_data(data_fn);
if ~isfield(data,'amplitudefl')
    errordlg('Data not found or not properly formatted','NIRFAST Error');
    error('Data not found or not properly formatted');
end
% remove zeroed data
ind = data.link(:,3)==0;
data.amplitudefl(ind,:) = []; clear ind
anom = log(data.amplitudefl);
% Only reconstructs fluorescence yield!

%*******************************************************
% load fine mesh for fwd solve: can input mesh structured variable
% or load from file
if ischar(fwd_mesh)==1
    fwd_mesh = load_mesh(fwd_mesh);
end
if ~strcmp(fwd_mesh.type,'fluor')
    errordlg('Mesh type is incorrect','NIRFAST Error');
    error('Mesh type is incorrect');
end
fwd_mesh.link = data.link;
clear data

etamuaf_sol=[output_fn '_etamuaf.sol'];

%**********************************************************
% Initiate log file

fid_log = fopen([output_fn '.log'],'w');
fprintf(fid_log,'Forward Mesh   = %s\n',fwd_mesh.name);
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
fprintf(fid_log,'Initial Guess muaf = %d\n',fwd_mesh.muaf(1));
% fprintf(fid_log,'Output Files   = %s',tau_sol);
fprintf(fid_log,'\n');

%***********************************************************
% get direct excitation field
% Flag mesh to not calculate the intrinsic emission and fluorescence
% emission fields
fwd_mesh.fl = 0; fwd_mesh.mm = 0; 
if isfield(fwd_mesh,'phix')~=0
    fwd_mesh = rmfield(fwd_mesh,'phix');
end
% calculate excitation field
data_fwd = femdata(fwd_mesh,frequency);
data_fwd.phi = data_fwd.phix;

%***********************************************************
% load recon_mesh
if ischar(recon_basis)
    recon_mesh = load_mesh(recon_basis);
    [fwd_mesh.fine2coarse,...
        recon_mesh.coarse2fine] = second_mesh_basis(fwd_mesh,recon_mesh);
elseif isstruct(recon_basis) == 0
    [fwd_mesh.fine2coarse,recon_mesh] = pixel_basis(recon_basis,fwd_mesh);
elseif isstruct(recon_basis) == 1
    if isfield(recon_basis,'nodes')
        recon_mesh = recon_basis;
        fwd_mesh.fine2coarse = recon_mesh.fine2coarse;
    else
        recon_mesh = recon_basis;
        [fwd_mesh.fine2coarse,...
            recon_mesh.coarse2fine] = second_mesh_basis(fwd_mesh,recon_mesh);
    end
end

%************************************************************
% initialize projection error
pj_error=[];

%*************************************************************
% modulation frequency
omega = 2*pi*frequency*1e6;
% set fluorescence variables
fwd_mesh.gamma = (fwd_mesh.eta.*fwd_mesh.muaf)./(1+(omega.*fwd_mesh.tau).^2);


% check for input regularization
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
%*************************************************************
% Calculate part of Jacobian which does not change at each iteration
% (call it "pre-Jacobian")
[Jpre,datafl,MASS_m] = prejacobian_fl(fwd_mesh,frequency,data_fwd);

%*************************************************************
% Iterate
for it = 1 : iteration
    
    % Update Jacobian with fluroescence field (changes at each iteration) 
    if it == 1
        [Jwholem,junk] = update_jacobian_fl(Jpre,fwd_mesh,frequency,data_fwd,MASS_m);
        clear junk
    else
        [Jwholem,datafl] = update_jacobian_fl(Jpre,fwd_mesh,frequency,data_fwd,MASS_m);
    end
    Jm = Jwholem.completem; clear Jwholem
    
    % Extract log amplitude reference data
    clear ref;
    ind = datafl.link(:,3)==0;
    datafl.amplitudem(ind,:) = []; clear ind
    ref(:,1) = log(datafl.amplitudem);
    
    % Calculate projection error
    data_diff = (anom-ref);
    pj_error = [pj_error sum(abs(data_diff.^2))];
        
    %***********************
    % Screen and Log Info
    disp('---------------------------------');
    disp(['Iteration_fl Number          = ' num2str(it)]);
    disp(['Projection_fl error          = ' num2str(pj_error(end))]);
    
    
    fprintf(fid_log,'---------------------------------\n');
    fprintf(fid_log,'Iteration_fl Number          = %d\n',it);
    fprintf(fid_log,'Projection_fl error          = %f\n',pj_error(end));
    
    if it ~= 1
        p = (pj_error(end-1)-pj_error(end))*100/pj_error(end-1);
        disp(['Projection error change   = ' num2str(p) '%']);
        fprintf(fid_log,'Projection error change   = %f %%\n',p);
        if (p) <= 2
            disp('---------------------------------');
            disp('STOPPING CRITERIA FOR FLUORESCENCE COMPONENT REACHED');
            fprintf(fid_log,'---------------------------------\n');
            fprintf(fid_log,'STOPPING CRITERIA FOR FLUORESCENCE COMPONENT REACHED\n');
            % set output
            data_recon.elements = fwd_mesh.elements;
            data_recon.etamuaf = fwd_mesh.etamuaf;
            break
        end
    end
    %*************************
    clear data_recon
       
    % Interpolate Jacobian onto recon mesh
    [Jm,recon_mesh] = interpolatef2r_fl(fwd_mesh,recon_mesh,Jm);
    Jm = Jm(:, 1:end/2); % take only intensity portion
    
    % Normalize Jacobian wrt fl source gamma
    Jm = Jm*diag([recon_mesh.gamma]);
    
    if strcmp(lambda.type, 'JJt')
        % build Hessian
        [nrow,ncol]=size(Jm);
        Hess = zeros(nrow);
        Hess = Jm*Jm';
       
        % add regularization
        reg = lambda.value.*(max(diag(Hess)));
        disp(['Regularization Fluor           = ' num2str(reg)]);
        fprintf(fid_log,'Regularization Fluor            = %f\n',reg);
        Hess = Hess+(eye(nrow).*reg);
        
        % Calculate update
        u = Jm'*(Hess\data_diff);
        u = u.*[recon_mesh.gamma];
    else
        % build Hessian
        [nrow,ncol]=size(Jm);
        Hess = zeros(ncol);
        Hess = Jm'*Jm;
        
        % add regularization
        reg = lambda.value.*(max(diag(Hess)));
        disp(['Regularization Fluor           = ' num2str(reg)]);
        fprintf(fid_log,'Regularization Fluor            = %f\n',reg);
        for i = 1 : ncol
            Hess(i,i) = Hess(i,i) + reg;
        end
        
        % Calculate update
        u = Hess\Jm'*data_diff;
        u = u.*[recon_mesh.gamma];
    end
    
    % value update:
    recon_mesh.gamma = recon_mesh.gamma+u;
    recon_mesh.etamuaf = recon_mesh.gamma.*(1+(omega.*recon_mesh.tau).^2);
    
    % negative constraint
    neg = find(recon_mesh.etamuaf <= 0);
    if isempty(neg) ~= 1
        recon_mesh.etamuaf(neg) = 10^-20;
    end
    
    % assuming we know eta
    recon_mesh.muaf = recon_mesh.etamuaf./recon_mesh.eta;
    clear u Hess Hess_norm tmp data_diff G
    
    % interpolate onto fine mesh
    [fwd_mesh,recon_mesh] = interpolatep2f_fl(fwd_mesh,recon_mesh);
    
    % filter
    if filter_n ~= 0
        disp('Filtering');
        fwd_mesh = mean_filter(fwd_mesh,filter_n);
    end
    
    %plotimage(fwd_mesh,fwd_mesh.eta.*fwd_mesh.muaf);
    %**********************************************************
    % Write solution to file
    
    if it == 1
        fid = fopen(etamuaf_sol,'w');
    else
        fid = fopen(etamuaf_sol,'a');
    end
    fprintf(fid,'solution %d ',it);
    fprintf(fid,'-size=%g ',length(fwd_mesh.nodes));
    fprintf(fid,'-components=1 ');
    fprintf(fid,'-type=nodal\n');
    fprintf(fid,'%g ',fwd_mesh.etamuaf);
    fprintf(fid,'\n');
    fclose(fid);
    
end
fin_it = it-1;
fclose(fid_log);
% Output recon basis mesh to use in subsequent reconstruction attempts.
recon_mesh.fine2coarse = fwd_mesh.fine2coarse;
fwd_mesh.recon_mesh = rmfield(recon_mesh,{'gamma','etamuaf','muaf','eta','tau'});


%******************************************************
% Sub functions
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
        %val_int(:,recon_mesh.elements(recon_mesh.coarse2fine(i,1),:)+NNC) = ...
        %    val_int(:,recon_mesh.elements(recon_mesh.coarse2fine(i,1),:)+NNC) + ...
        %    val(:,i+NNF)*recon_mesh.coarse2fine(i,2:end);
    elseif recon_mesh.coarse2fine(i,1) == 0
        dist = distance(fwd_mesh.nodes,fwd_mesh.bndvtx,recon_mesh.nodes(i,:));
        mindist = find(dist==min(dist));
        mindist = mindist(1);
        val_int(:,i) = val(:,mindist);
        %val_int(:,i+NNC) = val(:,mindist+NNF);
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
