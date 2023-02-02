function out = getspectra( lambda )
% % Hemoglobin
% http://omlc.ogi.edu/spectra/hemoglobin/
% scott prahl

% % Fat
% by R.L.P. van Veen and H.J.C.M. Sterenborg, A. Pifferi, A. Torricelli and R. Cubeddu
% http://omlc.ogi.edu/spectra/fat/

% % Water
% D. J. Segelstein, "The complex refractive index of water," University of
% Missouri-Kansas City, (1981).

    try
        load([fileparts(which('nirs.media.getspectra')) filesep 'spectra.mat'])
    catch
        load('NIRS_spectra.mat');
    end
    if(iscell(lambda))
        lambda=str2num(cell2mat(lambda));
    end
    
    
    out(:,1) = interp1(hbo(:,1),hbo(:,2),lambda);       % HbO extinction
    out(:,2) = interp1(hbr(:,1),hbr(:,2),lambda);       % HbR extinction
    out(:,3) = interp1(water(:,1),water(:,2),lambda); 	% water mua
    out(:,4) = interp1(fat(:,1),fat(:,2),lambda);    	% fat mua
    %out(:,5) = interp1(CytC(:,1),CytC(:,2),lambda);    	% CytC mua
    %out(:,6) = 519/10 * (lambda/500).^-3;             % MIE scattering
    

end