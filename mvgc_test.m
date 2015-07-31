%raw  = nirs.io.loadDirectory( '/Volumes/TinyDrive/Fetch/fnirs_data', {} );
raw  = nirs.io.loadDirectory( '~/Desktop/Tangos/data/processed', {} );

j = nirs.modules.RemoveStimless(); 
j = nirs.modules.Resample(j); 
j.Fs = 2;

j = nirs.modules.OpticalDensity(j); 
j = nirs.modules.BeerLambertLaw(j); 

hb = j.run(raw);

G1 = []; G2 = [];
for i = 3:100
    idx = randperm(length(hb),2);
    
    Y1 = hb(idx(1)).data(:,1:2:end);
    Y2 = hb(idx(2)).data(:,1:2:end);
    
    n = min(size(Y1,1), size(Y2,1));
    
    Y = [Y1(1:n,:) Y2(1:n,:)];
    
%     R(:,:,i) = corr(Y);
%     iR(:,:,i) = pinv(R(:,:,i));
    
    G1(:,:,i) = nirs.math.mvgc(Y, 16);
    
%     G2(:,:,i) = nirs.math.pwgranger(Y, 16);
    disp(i)
end

figure; imagesc(mean(G1,3)); colorbar; figure; imagesc(mean(G2,3)'), colorbar;