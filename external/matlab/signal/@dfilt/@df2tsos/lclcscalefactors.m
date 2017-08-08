function c = lclcscalefactors(Hd,c,nb,opts)
%LCLCSCALEFACTORS   Local cumulative scale factor computation.


%   Author(s): Dr. Guenter Dehner, R. Losada
%   Copyright 1988-2006 The MathWorks, Inc.

over = opts.over;
smax = opts.smax;
pnorm = opts.pnorm;


bb = 1;   aa = 1;

soc = Hd.sosMatrix;
for i= 1:nb,
    Hdpartial = cumsec(Hd,i);
    q1 = norm(Hdpartial,pnorm);

    switch over
        case {'wrap'}
            c(i) = smax/q1;
        case {'saturate'}
            Hdpartial = cumsec(Hd,i,true);

            Hdpartial.sosMatrix(1,1:2) = [(soc(i,4)*soc(i,2)-soc(i,5))...
                (soc(i,4)*soc(i,3)-soc(i,6))];
            q3 = norm(Hdpartial,pnorm);

            % Reuse temp object for speed
            Hdpartial.sosMatrix(1,1:2) = [(soc(i,4)*soc(i,3)-soc(i,6))...
                (soc(i,5)*soc(i,3)-soc(i,6)*soc(i,2))];
            q4 = norm(Hdpartial,pnorm);

            sp = [q1,q3,q4];
            maxsp = max(sp);
            c(i) = smax/maxsp;
        case {'satall'}
            a0 = aa;
            b0 = bb;
            aa = conv(aa,soc(i,4:6));
            bb = conv(bb,soc(i,1:3));
            b3 = conv(b0,[(soc(i,4)*soc(i,2)-soc(i,5))...
                (soc(i,4)*soc(i,3)-soc(i,6)) 0]);
            q3 = norm(dfilt.df2t(b3,aa),pnorm);
            b4 = conv(b0,[(soc(i,4)*soc(i,3)-soc(i,6))...
                (soc(i,5)*soc(i,3)-soc(i,6)*soc(i,2)) 0]);
            q4 = norm(dfilt.df2t(b4,aa),pnorm);
            b5 = conv(b0,[(soc(i,4)*soc(i,2))...
                (soc(i,4)*soc(i,3)+soc(i,5)*soc(i,2)-soc(i,6))...
                (soc(i,5)*soc(i,3))]);
            q5 = norm(dfilt.df2t(b5,aa),pnorm);
            q0 = norm(dfilt.df2t(b0,a0),pnorm);
            asoc= abs(soc(i,:));
            sp = [([1,asoc(4:6)])*q1,([1,asoc(2:3)])*q0,q3,q4,q5];
            maxsp = max(sp);
            c(1,i) = smax/maxsp;  
    end
end


% [EOF]
