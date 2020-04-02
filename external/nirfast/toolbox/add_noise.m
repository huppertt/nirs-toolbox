function data_noise = add_noise(data,amp,ph,fn)

% data_noise = add_noise(data,amp,ph,fn)
%
% Add normally distributed noise to data
%
% data is the filename of saved boundary data, or the workspace variable
% amp is the percentage noise in amplitude
% ph is the percentage noise in phase
% and fn is the name of the new file to be saved
 


% Initialize RANDN to a different state each time.
randn('state',sum(100*clock));

if ischar(data) ~= 0
    data = load_data(data);
end


if (isfield(data,'paa'))
    [nr,nc]=size(data.paa);
elseif isfield(data,'paaxflmm')
    [nr,nc]=size(data.paaxflmm);
elseif isfield(data,'paafl')
    [nr,nc]=size(data.paafl);
elseif isfield(data,'paax')
    [nr,nc]=size(data.paax);
elseif isfield(data,'paamm')
    [nr,nc]=size(data.paamm);
else
    disp('No appropriate data found.');
    return
end

% amplitude
n = randn(nr,nc/2);
%n = n./max(max(abs(n))); % make sure we go from -1 : 1
if (isfield(data,'paa'))
    data.paa(:,1:2:end) = data.paa(:,1:2:end) + ...
        amp.*(data.paa(:,1:2:end)./100).*n;
else
    if isfield(data,'paaxflmm')
        data.paaxflmm(:,1:2:end) = data.paaxflmm(:,1:2:end) + ...
            amp.*(data.paaxflmm(:,1:2:end)./100).*n;
        data.paax = data.paaxflmm(:,1:2);
        data.amplitudex = data.paax(:,1);
        data.paafl = data.paaxflmm(:,3:4);
        data.amplitudefl = data.paafl(:,1);
        data.paamm = data.paaxflmm(:,5:6);
        data.amplitudemm = data.paamm(:,1);
    else
        if isfield(data,'paafl')
            data.paafl(:,1:2:end) = data.paafl(:,1:2:end) + ...
                amp.*(data.paafl(:,1:2:end)./100).*n;
            data.amplitudefl = data.paafl(:,1);
        end
        if isfield(data,'paamm')
            data.paamm(:,1:2:end) = data.paamm(:,1:2:end) + ...
                amp.*(data.paamm(:,1:2:end)./100).*n;
            data.amplitudemm = data.paamm(:,1);
        end
        if isfield(data,'paax')
            data.paax(:,1:2:end) = data.paax(:,1:2:end) + ...
                amp.*(data.paax(:,1:2:end)./100).*n;
            data.amplitudex = data.paax(:,1);
        end
    end
end

%Phase
n = randn(nr,nc/2);
%n = n./max(max(abs(n))); % make sure we go from -1 : 1
if (isfield(data,'paa'))
    data.paa(:,2:2:end) = data.paa(:,2:2:end) + ph.*n;
else
    if isfield(data,'paaxflmm')
        data.paaxflmm(:,2:2:end) = data.paaxflmm(:,2:2:end) + ph.*n;
        data.paax = data.paaxflmm(:,1:2);
        data.phasex = data.paax(:,2);
        data.paafl = data.paaxflmm(:,3:4);
        data.phasefl = data.paafl(:,2);
        data.paamm = data.paaxflmm(:,5:6);
        data.phasemm = data.paamm(:,2);
    else
        if isfield(data,'paafl')
            data.paafl(:,2:2:end) = data.paafl(:,2:2:end) + ph.*n;
            data.phasefl = data.paafl(:,2);
        end
        if isfield(data,'paamm')
            data.paamm(:,2:2:end) = data.paamm(:,2:2:end) + ph.*n;
            data.phasemm = data.paamm(:,2);
        end
        if isfield(data,'paax')
            data.paax(:,1:2:end) = data.paax(:,1:2:end) + ph.*n;
            data.phasex = data.paax(:,2);
        end
    end
end

data_noise = data;
if nargin == 4
    save_data(data_noise,fn);
end