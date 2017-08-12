function odata = analyze_data(idata,mesh,tolamp,tolph,upper,lower)

% odata = analyze_data(idata,mesh,tolamp,tolph,upper,lower)
%
% Plots log amplitude & phase against distance (source to detector)
% to find bad data points and replace with NaNs
%
% idata is the structured data variable
% odata is the output, with bad points removed
% mesh is the structured mesh variable
% tolamp is the tolerance off the line of fit (amplitude) for data removal
%       can be a vector for spectral
% tolph is the tolerance off the line of fit (phase) for data removal
%       can be a vector for spectral
% upper is an upper limit for source/detector distance
%       can be a vector for spectral
% lower is a lower limit for source/detector distance
%       can be a vector for spectral


% error checking
if tolamp < 0
    errordlg('Amplitude tolerance must be nonnegative','NIRFAST Error');
    error('Amplitude tolerance must be nonnegative');
end
if tolph < 0
    errordlg('Phase tolerance must be nonnegative','NIRFAST Error');
    error('Phase tolerance must be nonnegative');
end
if upper < lower
    errordlg('Upper limit on source/detector distance must be greater than lower limit','NIRFAST Error');
    error('Upper limit on source/detector distance must be greater than lower limit');
end

odata = idata;

% load data
if ischar(idata) ~= 0
    idata = load_data(idata);
end

% load mesh
if ischar(mesh) ~= 0
    mesh = load_mesh(mesh);
end

% find data
phasedata = 0;
if (strcmp(mesh.type,'stnd') || strcmp(mesh.type,'spec') || strcmp(mesh.type,'stnd_bem') || strcmp(mesh.type,'spec_bem')) ...
        && isfield(idata,'paa')
    data_big = idata.paa;
    phasedata = 1;
elseif (strcmp(mesh.type,'fluor') || strcmp(mesh.type,'fluor_bem')) ...
        && (isfield(idata,'paafl') || isfield(idata,'amplitudefl'))
    if isfield(idata,'paafl')
        data_big = idata.paafl;
        phasedata = 1;
    else
        data_big = idata.amplitudefl;
    end
else
    errordlg('There is no data or it is not properly formatted','NIRFAST Error');
    error('There is no data or it is not properly formatted');
end

% source/detector distances
dist = zeros(length(mesh.link),1);
for i = 1:length(mesh.link)
    snum = mesh.link(i,1);
    mnum = mesh.link(i,2);
    snum = mesh.source.num == snum;
    mnum = mesh.meas.num == mnum;
    if sum(snum)==0 || sum(mnum)==0
        dist_full(i,1)=0;
        mesh.link(i,3)=0;
    else
        dist_full(i,1) = sqrt(sum((mesh.source.coord(snum,:) - ...
        mesh.meas.coord(mnum,:)).^2,2)); 
    end
end

% get an index from link file of data to actually use
linki = logical(mesh.link(:,3));


for i=1:1+phasedata:size(data_big,2)
    
    if ~isscalar(upper) && (i+1)/2 > length(upper) || ...
            ~isscalar(lower) && (i+1)/2 > length(lower) || ...
            ~isscalar(upper) && (i+1)/2 > length(upper) || ...
            ~isscalar(upper) && (i+1)/2 > length(upper)
        errordlg('The mesh and data wavelength information do not match','NIRFAST Error');
        error('The mesh and data wavelength information do not match');
    end
    
    if phasedata
        data = data_big(:,i:i+1);
    else
        data = data_big(:,i);
    end
    
    % Set lnrI, ph
    data_tmp = data(linki,:);
    dist = dist_full(linki);
    [j,k] = size(data_tmp(:,1));
    [j2,k2] = size(dist);
    lnrI = log(data_tmp(:,1).*dist);
    if phasedata
        ph = data_tmp(:,2);
    end
    
    % remove data based on upper and lower bounds (just for the fitting
    % here)
    lnrI_tmp = lnrI;
    if phasedata
        ph_tmp = ph;
    end
    if isscalar(upper)
        up = upper;
    else
        up = upper((i+1)/2);
    end
    if isscalar(lower)
        lo = lower;
    else
        lo = lower((i+1)/2);
    end
    
    for k=1:1:j
        if dist(k) > up || dist(k) < lo
            lnrI_tmp(k)=NaN;
            if phasedata
                ph_tmp(k)=NaN;
            end
        end   
    end

    % Calculate best fit line
    m1 = polyfit(dist(isnan(lnrI_tmp)==0),lnrI(isnan(lnrI_tmp)==0),1);
    if phasedata
        m2 = polyfit(dist(isnan(ph_tmp)==0),ph(isnan(ph_tmp)==0),1);
    end
    
    lnrI = log(data(:,1).*dist_full);
    if phasedata
        ph = data(:,2);
    end

    % Remove points too far from the line of fit
    badpoints_amp = zeros(length(dist_full),1);
    badpoints_ph = zeros(length(dist_full),1);
    if isscalar(tolamp)
        ta = tolamp;
    else
        ta = tolamp((i+1)/2);
    end
    if isscalar(tolph)
        tp = tolph;
    else
        tp = tolph((i+1)/2);
    end
    
    for k=1:1:length(dist_full)
        if abs((lnrI(k) - (m1(1)*dist_full(k)+m1(2)))) > ta...
                || dist_full(k) > up || dist_full(k) < lo
            badpoints_amp(k)=1;
        end
        if phasedata &&...
                (abs((ph(k) - (m2(1)*dist_full(k)+m2(2)))) > tp...
                || dist_full(k) > up || dist_full(k) < lo )
            badpoints_ph(k)=1;
        end     
    end
    
    % set NaNs
    if phasedata
        odata.link(badpoints_amp==1,2+(i+1)/2) = 0;
        odata.link(badpoints_ph==1,2+(i+1)/2) = 0;
    else
        odata.link(badpoints_amp==1,2+i) = 0;
    end
end
    