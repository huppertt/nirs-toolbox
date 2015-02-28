function data = load_data(fn, wv_array)

% data = load_data(fn, wv_array)
%
% Loads data from the file in the appropriate data format
%
% fn is the filename
% wv_array is optional - wavelength array for spectral meshes
%       if wv_array is specified, data will be formatted for spectral recon
% data is the data read
%
% spectral info should have the wavelengths on the first line, like:
%
% w661 w808 w910
%
% fluorescence should have whichever of the following labels are relevant
% at the top:
%
% xphase   xamplitude   flphase   flamplitude   mphase   mamplitude

data.name = 'nirfast';

if ischar(fn) ~= 0
    
    if ~exist(fn,'file')
        errordlg('The data file could not be found','NIRFAST Error');
        error('The data file could not be found');
    end
    
    datatemp = importdata(fn);
    
    % support for new matlab importdata format
    if isfield(datatemp,'textdata') && size(datatemp.textdata,2) > 1
        tempvar = [];
        for wv=1:1:size(datatemp.textdata,2)
            tempvar = [tempvar ' ' char(datatemp.textdata(wv))];
        end
        datatemp.textdata = {tempvar(2:end)};
    end
    
    if isfield(datatemp,'textdata')
        text = datatemp.textdata;
        test = char(text(1));
        datatemp = datatemp.data;
    else
        test = ' ';
    end
    
    if strcmp(test(1),'s') == 1
        if isempty(strfind(test,'active')) == 0
            test(1:strfind(test,'active')+5)=[];
            data.link = datatemp(:,1:3);
            datatemp = datatemp(:,4:end);
        elseif isempty(strfind(test,'w')) == 0
            foo = strfind(test,'w');
            test(1:foo(1)-1) = [];
            data.link = datatemp(:,1:2);
            datatemp = datatemp(:,3:end);
        end
        test = test(find(isspace(test)==0):end);
    end
    
    % FLUORESCENCE
    if (strcmp(test(1),'x') == 1 || strcmp(test(1),'f') == 1 || strcmp(test(1),'m') == 1)
        
        % find labels
        S = test;
        i = 1;
        while ~isempty(S)
            [T,S] = strtok(S);
            if strcmp(T,'xphase')
                data.phasex = datatemp(:,i);
            elseif strcmp(T,'xamplitude')
                data.amplitudex = datatemp(:,i);
            elseif strcmp(T,'mmphase')
                data.phasemm = datatemp(:,i);
            elseif strcmp(T,'mmamplitude')
                data.amplitudemm = datatemp(:,i);
            elseif strcmp(T,'flphase')
                data.phasefl = datatemp(:,i);
            elseif strcmp(T,'flamplitude')
                data.amplitudefl = datatemp(:,i);
            end
            i = i + 1;
        end
        
        % construct paa*
        if isfield(data,'phasex') && isfield(data,'amplitudex')
            data.paax = [data.amplitudex data.phasex];
        end
        if isfield(data,'phasemm') && isfield(data,'amplitudemm')
            data.paamm = [data.amplitudemm data.phasemm];
        end
        if isfield(data,'phasefl') && isfield(data,'amplitudefl')
            data.paafl = [data.amplitudefl data.phasefl];
        end
        if isfield(data,'paax') && isfield(data,'paafl')
            data.paaxfl = [data.paax data.paafl];
        end
        if isfield(data,'paaxfl') && isfield(data,'paamm')
            data.paaxflmm = [data.paaxfl data.paamm];
        end
        
        
        % SPECTRAL
    elseif isempty(strfind(test,'w')) == 0
        S = char(text);
        wloc = strfind(S,'w');
        for i=1:1:numel(wloc)-1
            data.wv(i) = str2num(strtrim(S(wloc(i)+1:wloc(i+1)-1)));
        end
        data.wv(end+1) = str2num(strtrim(S(wloc(end)+1:end)));
        
        % finish writing link
        data.link = [data.link datatemp(:,1:3:end)];
        datatemp(:,1:3:end) = [];
        data.paa = datatemp;
              
        
        % STANDARD
    else
        [junk,n] = size(datatemp);
        data.paa = datatemp;
        if n>1
            data.phase = data.paa(:,2);
            data.amplitude = data.paa(:,1);
        else
            data.amplitude = data.paa(:,1);
        end
    end
    
else
    data = fn;
end

% if wv_array exists, only take those wavelengths
% modified for consistency, so that it still outputs Amplitude and phase
% in radians, and only the wavelength at which it is. HD 07-09-09
% fixed as data.paa for even wavelengths were wrong HD 07-09-09
if (nargin > 1 && isfield(data,'wv') && isfield(data,'link'))
    link = data.link(:,1:2);
    for i = 1:length(wv_array)
        link = [link data.link(:,find(data.wv == wv_array(i))+2)];
    end
    data.link = link;
end

if (nargin > 1 && isfield(data,'wv'))
    anom_big = [];
    for i = 1:length(wv_array)
        index_wv = find(data.wv == wv_array(i));
        if isempty(index_wv)
            display([ num2str(wv_array(i)) ' nm data not available']);
            data = [];
            return;
        end
        anom(:,1) = data.paa(:,(index_wv(1)-1)*2+1);
        anom(:,2) = data.paa(:,(index_wv(1)-1)*2+2);
        %anom(:,1) = log(anom(:,1));
        %anom(:,2) = anom(:,2)/180.0*pi;
        %anom(find(anom(:,2)<0),2) = anom(find(anom(:,2)<0),2) + (2*pi);
        %anom(find(anom(:,2)>(2*pi)),2) = anom(find(anom(:,2)>(2*pi)),2) - (2*pi);
        %anom = reshape(anom',length(anom)*2,1);
        anom_big = [anom_big anom];
        clear anom
    end
    data.paa = anom_big;
    data.wv = wv_array;
end