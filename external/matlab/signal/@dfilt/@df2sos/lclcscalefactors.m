function c = lclcscalefactors(Hd,c,nb,opts)
%LCLCSCALEFACTORS   Local cumulative scale factor computation.


%   Author(s): Dr. Guenter Dehner, R. Losada
%   Copyright 2003-2006 The MathWorks, Inc.

over = opts.over;
smax = opts.smax;
pnorm = opts.pnorm;
svm = opts.ScaleValueConstraint;

bb = 1;   aa = 1;
soc = Hd.sosMatrix;

for i= 1:nb
    Hdpartial = cumsec(Hd,i,true);
    q1 = norm(Hdpartial,pnorm);
    switch over
        case {'wrap' 'saturate'}
            c(1,i) = smax/q1;
            
            if ~strcmpi(svm,'unit'),
                Hdpartial = cumsec(Hd,i,false);
                c(2,i) = smax/norm(Hdpartial,pnorm);
            end
        case {'satall'}
            aa = conv(aa,soc(i,4:6));
            b2 = conv(bb,soc(i,1:3));
            b3 = conv(bb,[0 -soc(i,5) -soc(i,6)]);
            htemp = dfilt.df2(b3,aa);
            q3 = norm(htemp,pnorm);
            sp1 =[q1,abs(soc(i,5))*q1,abs(soc(i,6))*q1,q3,soc(i,4)*q1];
            
            htemp = dfilt.df2(b2,aa);
            q2 = norm(htemp,pnorm);
            b4 = conv(bb,[0 soc(i,2) soc(i,3)]);
            htemp = dfilt.df2(b4,aa);
            q4 = norm(htemp,pnorm);
            sp2 = [q2,q1,abs(soc(i,2))*q1,abs(soc(i,3))*q1,q4];
            
            if strcmpi(svm,'unit'),
                sp = [sp1 sp2];  
                maxsp = max (sp);
                c(1,i) = smax/maxsp;
            else
                maxsp1 = max(sp1);
                c(1,i) = smax/maxsp1; 
                maxsp2 = max(sp2);
                c(2,i) = smax/maxsp2; 
            end
            bb = b2;
    end
end
if strcmpi(svm,'unit'),
    Hdpartial = cumsec(Hd,nb,false);
    c(2,nb) = smax/norm(Hdpartial,pnorm);    
end


% [EOF]
