function c = lclcscalefactors(Hd,c,nb,opts)
%LCLCSCALEFACTORS   Local cumulative scale factor computation.


%   Author(s): Dr. Guenter Dehner, R. Losada
%   Copyright 1988-2006 The MathWorks, Inc.

over = opts.over;
smax = opts.smax;
pnorm = opts.pnorm;

bb = 1;   aa = 1;
soc = Hd.sosMatrix;
q0 = 1;
for i= 1:nb,
    Hdpartial = cumsec(Hd,i);
    q1 = norm(Hdpartial,pnorm);
    switch over
        case {'wrap' 'saturate'}
            c(1,i) = smax/q1;
        case {'satall'}  % all multiplier and adder outputs are checked
            a0 = aa;
            b0 = bb;
            aa = conv(aa,soc(i,4:6));
            bb = conv(bb,soc(i,1:3));
            sp = [1,soc(i,4),abs(soc(i,5)),abs(soc(i,6))]*q1;
            sp = [sp,[1,abs(soc(i,2)),abs(soc(i,3))]*q0];
            b3 = conv(bb,[0 -soc(i,5) -soc(i,6)]);
            b4 = conv(b0,[0 soc(i,2) soc(i,3)]);
            b5 = conv(b0,[1 soc(i,2) soc(i,3)]);
            sp = [sp,norm(dfilt.df1(b3,aa),pnorm),norm(dfilt.df1(b4,a0),pnorm),...
                norm(dfilt.df1(b5,a0),pnorm)];
            maxsp = max(sp);
            c(1,i) = smax/maxsp;
    end
    q0 = q1;   % output = input of next block
end


% [EOF]
