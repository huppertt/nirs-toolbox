function [fwd_mesh,pj_error] = reconstruct_fl_region(fwd_mesh,...
    frequency,...
    data_fn,...
    iteration,...
    lambda,...
    output_fn,...
    filter_n,...
    region)

% [fwd_mesh,pj_error] = reconstruct_fl_region(fwd_mesh,...
%                                      frequency,...
%                                      data_fn,...
%                                      iteration,...
%                                      lambda,...
%                                      output_fn,...
%                                      filter_n,...
%                                      region)
%
% region-based reconstruction program for fluorescence meshes
%
% fwd_mesh is the input mesh (variable or filename)
% frequency is the modulation frequency (MHz)
% data_fn is the boundary data (variable or filename)
% iteration is the max number of iterations
% lambda is the initial regularization value
%    e.g. [20 20 15; 20 20 15] for a 3D mesh with 2 regions
% output_fn is the root output filename
% filter_n is the number of mean filters
% region is an array of the regions (e.g. [0 1 2])



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

%**********************************************************
% initialize projection error
pj_error=[];

%*************************************************************
% modulation frequency
omega = 2*pi*frequency*1e6;
% set fluorescence variables
fwd_mesh.gamma = (fwd_mesh.eta.*fwd_mesh.muaf)./(1+(omega.*fwd_mesh.tau).^2);

%*************************************************************
% Calculate region mapper
disp('calculating regions');
if ~exist('region','var') || (exist('region','var') && isempty(region))
    region = unique(fwd_mesh.region);
end
K = region_mapper(fwd_mesh,region);

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
    
    % Read reference data
    clear ref;
    ind = datafl.link(:,3)==0;
    datafl.amplitudem(ind,:) = []; clear ind
    ref(:,1) = log(datafl.amplitudem);
    
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
    
    Jm = Jm*diag([fwd_mesh.gamma]);
    
    % reduce J into regions!
    Jm = Jm*K;
    
    % build Hessian
    [nrow,ncol]=size(Jm);
    Hess = zeros(nrow);
    Hess = Jm*Jm';
    
    % add regularization
    reg = lambda.*(max(diag(Hess)));
    disp(['Regularization Fluor           = ' num2str(reg)]);
    fprintf(fid_log,'Regularization Fluor            = %f\n',reg);
    Hess = Hess+(eye(nrow).*reg);
    
    % Calculate update
    u = Jm'*(Hess\data_diff);
    
    % use region mapper to unregionize!
    u = K*u;
    u = u.*[fwd_mesh.gamma];
    
    % value update:
    fwd_mesh.gamma = fwd_mesh.gamma+u;
    fwd_mesh.etamuaf = fwd_mesh.gamma.*(1+(omega.*fwd_mesh.tau).^2);
    
    % negative constraint
    fwd_mesh.etamuaf(find(fwd_mesh.etamuaf <= 0)) = 10e-20;
    
    % assuming we know eta
    fwd_mesh.muaf = fwd_mesh.etamuaf./fwd_mesh.eta;
    clear u Hess Hess_norm tmp data_diff G
    
    % filter
    if filter_n ~= 0
        disp('Filtering');
        fwd_mesh = mean_filter(fwd_mesh,filter_n);
    end
    
    
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
