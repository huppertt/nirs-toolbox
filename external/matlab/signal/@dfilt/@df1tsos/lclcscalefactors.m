function c = lclcscalefactors(Hd,c,nb,opts)
%LCLCSCALEFACTORS   Local cumulative scale factor computation.


%   Author(s): Dr. Guenter Dehner, R. Losada
%   Copyright 2003-2006 The MathWorks, Inc.

over = opts.over;
smax = opts.smax;
pnorm = opts.pnorm;
svm = opts.ScaleValueConstraint;

soc = Hd.sosMatrix;
c=zeros(2,nb+1);

bb = 1;   aa = 1;

for i= 1:nb,

    Hdpartial = cumsec(Hd,i,false);
    q2 = norm(Hdpartial,pnorm);

    Hdpartial = cumsec(Hd,i,true); % Make sure to do this last, since Hdpartial is used later on
    q1 = norm(Hdpartial,pnorm);


    switch over
        case {'wrap' }
            c(1,i) = smax/q1;  
            if ~strcmpi(svm,'unit') || i==nb,
                c(2,i) = smax/q2; 
            end
        case {'saturate'}
            Hdpartial.sosMatrix(1,1:2) = [-soc(i,5) -soc(i,6)];
            q3 = norm(Hdpartial,pnorm);
            Hdpartial.sosMatrix(1,1:2) = [soc(i,2) soc(i,3)];
            q4 = norm(Hdpartial,pnorm);
            q5 = abs(soc(i,3))*q1;
            sp1 = [q1,q3];    
            maxsp1 = max(sp1);
            c(1,i) = smax/maxsp1; 
            if strcmpi(svm,'unit'),
                sp2 = [q4,q5]; 
            else
                sp2 = [q2,q4,q5];
            end
            maxsp2 = max(sp2);
            c(2,i) = smax/maxsp2; 
        case {'satall'}
            aa = conv(aa,soc(i,4:6));
            b2 = conv(bb,soc(i,1:3));
            b3 = conv(bb,[-soc(i,5) -soc(i,6) 0 ]);
            q3 = norm(dfilt.df1t(b3,aa),pnorm);
            sp1 = [q1,abs(soc(i,5))*q1,abs(soc(i,6))*q1,q3,soc(i,4)*q1];
            maxsp1 = max(sp1);
            c(1,i) = smax/maxsp1; 
            q2 = norm(dfilt.df1t(b2,aa),pnorm);
            q0 = max([1,abs(soc(i,2)),abs(soc(i,3))])*q1;
            b4 = conv(bb,[soc(i,2) soc(i,3) 0 ]);
            q4 = norm(dfilt.df1t(b4,aa),pnorm);
            sp2 = [q2,q1,abs(soc(i,2))*q1,abs(soc(i,3))*q1,q4];
            maxsp2 = max(sp2);
            c(2,i) = smax/maxsp2;
            bb = b2;
    end
end
if strcmpi(svm,'unit') && ~strcmpi(over,'wrap')
    for i=2:nb,
        c(1,i) = min([c(1,i),c(2,i-1)]);
        c(2,i-1)=0;        
    end
end

% [EOF]
