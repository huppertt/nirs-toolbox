function plot_data(data)

% plot_data(data)
%
% Plots phase and amplitude of the data
%
% data is the structured data variable


form = 'b.-';

% load from file if path given
if ischar(data)
    data = load_data(data);
end

% TIME RESOLVED
if isfield(data,'tpsf')
    figure
    imagesc(abs(data.tpsf));
    xlabel('time point');
    ylabel('source/detector pair');
    title('Time Resolved Data');
elseif isfield(data,'tpsfx')
    figure
    imagesc(abs(data.tpsfx));
    xlabel('time point');
    ylabel('source/detector pair');
    title('Time Resolved Data');
% STANDARD OR SPECTRAL
elseif isfield(data,'paa')
    [j,k] = size(data.paa);
    m = 1;
    for i=1:2:k
        figure;
        setfont(14);
        
        linki = logical(data.link(:,m+2));
        m = m + 1;
        semilogy(data.paa(linki,i),form);
        if isfield(data,'wv')
            title([num2str(data.wv(m-1)) 'nm Amplitude']);
        else
            title('Amplitude');
        end
        ylabel('Amplitude (W/mm^2)');
        xlabel('Measurement');
        
        figure
        setfont(14);
        
        plot(data.paa(linki,i+1),form);
        if isfield(data,'wv')
            title([num2str(data.wv(m-1)) 'nm Phase']);
        else
            title('Phase');
        end
        ylabel('Phase (deg)');
        xlabel('Measurement');
    end
% FLUORESCENCE
elseif isfield(data,'paaxflmm')
    linki = logical(data.link(:,3));
    [j,k] = size(data.paaxflmm);
    figure;
    setfont(14);
    
    semilogy(data.paax(linki,1),form);
    title('Excitation Amplitude');
    ylabel('Amplitude (W/mm^2)');
    xlabel('Measurement');
    figure;
    setfont(14);
    
    plot(data.paax(linki,2),form);
    title('Excitation Phase');
    ylabel('Phase (deg)');
    xlabel('Measurement');
    figure;
    setfont(14);
    
    semilogy(data.paafl(linki,1),form);
    title('Fluorescence Amplitude');
    ylabel('Amplitude (W/mm^2)');
    xlabel('Measurement');
    figure;
    setfont(14);
    
    plot(data.paafl(linki,2),form);
    title('Fluorescence Phase');
    ylabel('Phase (deg)');
    xlabel('Measurement');
    figure;
    setfont(14);
    
    semilogy(data.paamm(linki,1),form);
    title('Emission Amplitude');
    ylabel('Amplitude (W/mm^2)');
    xlabel('Measurement');
    figure;
    setfont(14);
    
    plot(data.paamm(linki,2),form);
    title('Emission Phase');
    ylabel('Phase (deg)');
    xlabel('Measurement');
elseif isfield(data,'paaxfl')
    linki = logical(data.link(:,3));
    [j,k] = size(data.paaxfl);
    figure;
    setfont(14);
    
    semilogy(data.paax(linki,1),form);
    title('Excitation Amplitude');
    ylabel('Amplitude (W/mm^2)');
    xlabel('Measurement');
    figure;
    setfont(14);
    
    plot(data.paax(linki,2),form);
    title('Excitation Phase');
    ylabel('Phase (deg)');
    xlabel('Measurement');
    figure;
    setfont(14);
    
    semilogy(data.paafl(linki,1),form);
    title('Fluorescence Amplitude');
    ylabel('Amplitude (W/mm^2)');
    xlabel('Measurement');
    figure;
    setfont(14);
    
    plot(data.paafl(linki,2),form);
    title('Fluorescence Phase');
    ylabel('Phase (deg)');
    xlabel('Measurement');
elseif isfield(data,'paax')
    linki = logical(data.link(:,3));
    [j,k] = size(data.paax);
    figure;
    setfont(14);
    
    semilogy(data.paax(linki,1),form);
    title('Excitation Amplitude');
    ylabel('Amplitude (W/mm^2)');
    xlabel('Measurement');
    figure;
    setfont(14);
    
    plot(data.paax(linki,2),form);
    title('Excitation Phase');
    ylabel('Phase (deg)');
    xlabel('Measurement');
elseif isfield(data,'paafl')
    linki = logical(data.link(:,3));
    [j,k] = size(data.paafl);
    figure;
    setfont(14);
    
    semilogy(data.paafl(linki,1),form);
    title('Fluorescence Amplitude');
    ylabel('Amplitude (W/mm^2)');
    xlabel('Measurement');
    figure;
    setfont(14);
    
    plot(data.paafl(linki,2),form);
    title('Fluorescence Phase');
    ylabel('Phase (deg)');
    xlabel('Measurement');
elseif isfield(data,'paamm')
    linki = logical(data.link(:,3));
    [j,k] = size(data.paamm);
    figure;
    setfont(14);
    
    semilogy(data.paamm(linki,1),form);
    title('Emission Amplitude');
    ylabel('Amplitude (W/mm^2)');
    xlabel('Measurement');
    figure;
    setfont(14);
    
    plot(data.paamm(linki,2),form);
    title('Emission Phase');
    ylabel('Phase (deg)');
    xlabel('Measurement');
elseif isfield(data,'amplitudefl')
    linki = logical(data.link(:,3));
    [j,k] = size(data.amplitudefl);
    figure;
    setfont(14);
    
    semilogy(data.amplitudefl(linki),form);
    title('Fluorescence Amplitude');
    ylabel('Amplitude (W/mm^2)');
    xlabel('Measurement');
else
    errordlg('Data not found or not properly formatted','NIRFAST Error');
    error('Data not found or not properly formatted');
end

end


function setfont(s)

title([],'fontsize',s);
xlabel([],'fontsize',s);
xlabel([],'fontsize',s);
set(gca,'fontsize',s);
        
end