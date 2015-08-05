function createFakeDataset( folder )
    % SD variable
    SD.SrcPos = [];
    SD.SrcPos(:,1) = [-80:20:-20 20:20:80];
    SD.SrcPos(:,2) = 20;
    SD.SrcPos(:,3) = 0;

    SD.DetPos = [];
    SD.DetPos(:,1) = [-70:20:-30 30:20:70];
    SD.DetPos(:,2) = 0;
    SD.DetPos(:,3) = 0;

    SD.NumSrc = size(SD.SrcPos,1);
    SD.NumDet = size(SD.DetPos,1);

    n = SD.NumSrc;
    SD.MeasList = [
        1 1;
        2 1;
        2 2;
        3 2;
        3 3;
        4 4;
        5 4;
        6 4;
        6 5;
        7 5;
        7 6;
        8 6];

    SD.MeasList(:,3:4) = 1;
    SD.MeasList = [SD.MeasList; SD.MeasList];

    SD.MeasList(end/2:end,4) = 2;
    
    % time variable
    t = (0:5*60*4-1)' / 4;
    
    % StimDesign
    on = 10:20:300-20;
    dur = 10*ones(size(on));
    amp = 10*ones(size(on));
    
    StimDesign(1).cond  = 'A';
    StimDesign(1).onset = on(1:2:end);
    StimDesign(1).dur   = dur(1:2:end);
    StimDesign(1).amp   = amp(1:2:end);
    
    StimDesign(2).cond  = 'B';
    StimDesign(2).onset = on(2:2:end);
    StimDesign(2).dur   = dur(2:2:end);
    StimDesign(2).amp   = amp(2:2:end);
end

function Y = simulateNoise(ntime, nchan)
    Y = randn(ntime, nchan);
    
    f = [1 -0.5 -0.3 0 0 0.3 0 0 -0.1];
    Y = filter(1, f, Y);
end