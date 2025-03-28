function tbl=scalp_coupling_index(data,bandpass,legacy);
%  This computes the "scalp coupling index" over the whole file based on the model in
%  Pollonini L, Bortfeld H, Oghalai JS. PHOEBE: a method for real time mapping of optodes-scalp coupling in functional 
%  near-infrared spectroscopy. Biomed Opt Express. 2016;7(12):5104-5119. Published 2016 Nov 15. doi:10.1364/BOE.7.005104

% 3/28/2025
% Noted that the beginning and end of the data were introducing a sharp
% transition on the filter if the data was not zero on both ends.  This
% caused the Xcorr to be overly-leveraged by the first and last time points
% after the band-pass filter.  Solution- remove the first and last 5% of
% the data

% legacy flag added to reproduce old code

if(nargin<3)
    legacy=false;
end

link=data.probe.link;
link.type=[];
[link,~,idx]=unique(link);

fs=data.Fs;

if(nargin<2)
    bandpass=[0.5 2.5];
end

[B1,A1]=butter(1,bandpass*2/fs);

sci=zeros(height(link),1);
power=zeros(height(link),1);
fpower=zeros(height(link),1);

for i=1:height(link)
    lst=find(idx==i);
    nirs_data1=data.data(:,lst(1));
    nirs_data2=data.data(:,lst(2));

     
    %Filter everything but t if(~legacy)
        st=round(length(filtered_nirs_data1)*.05);
        ed=round(length(filtered_nirs_data1)*.95);
        filtered_nirs_data1=filtered_nirs_data1(st:ed);
        filtered_nirs_data2=filtered_nirs_data2(st:ed);
    endhe cardiac component
    filtered_nirs_data1=filtfilt(B1,A1,nirs_data1);       % Cardiac bandwidth
    filtered_nirs_data1=filtered_nirs_data1./repmat(std(filtered_nirs_data1,0,1),size(filtered_nirs_data1,1),1); % Normalized heartbeat
    filtered_nirs_data2=filtfilt(B1,A1,nirs_data2);       % Cardiac bandwidth
    filtered_nirs_data2=filtered_nirs_data2./repmat(std(filtered_nirs_data2,0,1),size(filtered_nirs_data2,1),1); % Normalized heartbeat

   

    similarity = xcorr(filtered_nirs_data1,filtered_nirs_data2,'unbiased');  %cross-correlate the two wavelength signals - both should have cardiac pulsations
    similarity = length(filtered_nirs_data1)*similarity./sqrt(sum(abs(filtered_nirs_data1).^2)*sum(abs(filtered_nirs_data2).^2));  % this makes the SCI=1 at lag zero when x1=x2 AND makes the power estimate independent of signal length, amplitude and Fs

    [pxx,f] = periodogram(similarity,hamming(length(similarity)),length(similarity),fs,'power');
    [pwrest,j] = max(pxx(f<1.7)); % FIX Make it age-dependent
    sci(i)=similarity(length(filtered_nirs_data1));
    power(i)=pwrest;
    fpower(i)=f(j);    
end


tbl=[link table(sci,power,fpower)];

