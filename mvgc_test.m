raw  = nirs.io.loadDirectory( '.', {} );

j = nirs.modules.RemoveStimless(); 
j = nirs.modules.Resample(j); 
j.Fs = 4;
j = nirs.modules.OpticalDensity(j); 
j = nirs.modules.BeerLambertLaw(j); 

hb = j.run(raw);

G = []; R = []; iR = [];

for i = 3:100
    idx = randperm(length(hb),2);
    
    Y1 = hb(idx(1)).data(:,1:2:end);
    Y2 = hb(idx(2)).data(:,1:2:end);
    
    n = min(size(Y1,1), size(Y2,1));
    
    Y = [Y1(1:n,:) Y2(1:n,:)];
    
%     R(:,:,i) = corr(Y);
%     iR(:,:,i) = pinv(R(:,:,i));
    
    G(:,:,i) = nirs.math.mvgc(Y, 16);
%     V(:,:,i) = cov(Y);
end