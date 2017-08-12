function plot_emiss_spec(mesh,data,true_fl_fn,source,meas)

% plot_emiss_spec(mesh,data,true_fl_fn,source,meas)
%
% plots the emission spectrum from spectral emission forward data
%
% mesh is the input mesh (variable or filename)
% data is the spectral emission data (variable or filename)
% true_fl_fn is the drug file
% source is the source number
% meas is the detector number

if ~exist(true_fl_fn,'file')
    errordlg('Drug file not found','NIRFAST Error');
    error('Drug file not found');
end

if ischar(mesh) == 1
    mesh = load_mesh(mesh);
end
if ischar(data) == 1
    data = load_data(data);
end
if ~isfield(data,'paa') || ~isfield(data,'wv')
    errordlg('Data not found or not properly formatted','NIRFAST Error');
    error('Data not found or not properly formatted');
end

figure
index = find(data.link(:,1)==source & data.link(:,2)==meas);

true_fl = load(true_fl_fn);
A = plot(true_fl(:,1),true_fl(:,2)/max(true_fl(:,2)),data.wv,...
    data.paa(index,1:2:end-1)/max(data.paa(index,1:2:end-1)));

ylabel('Normalized Intensity','FontSize',16);
xlabel('Wavelength (nm)','FontSize',16);
legend('Fluor in Dilute Solution','Fluor through tissue');

set(A(1),'LineStyle',':','LineWidth',3,'Color','r');
set(A(2),'Marker','d','LineWidth',3);