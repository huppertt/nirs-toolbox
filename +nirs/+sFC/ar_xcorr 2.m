function [R,p,dfe]=ar_xcorr(data,modelorder,maxlags,robust_flag,avg_method)
if(nargin<5)
    avg_method = 'mean';
end

if(nargin<4)
    robust_flag=true;
end

if(nargin<3)
    maxlags=[];
end

if(nargin<2 || isempty(modelorder))
    modelorder=20;
end

if(~isempty(strfind(class(data),'.core.Data')))
    Fs=data.Fs;
    data=data.data;
else
    Fs=1;
end

if(isstr(modelorder))
    p = Fs*str2num(modelorder(1:strfind(modelorder,'x')-1));
else
    p=modelorder;
end

if(isstr(maxlags))
    maxlags = ceil(Fs*str2num(maxlags(1:strfind(maxlags,'x')-1)));
end

if(~isreal(data))
    mask=(imag(data)>0);
    
else
    mask=ones(size(data));
end

[yfilt,f] = nirs.math.innovations(real(data),p);

% Mask out boundary values
for ch = 1:size(yfilt,2)
    yfilt(1:length(f{ch})) = nan;
end

if(robust_flag)
    [R,p]=nirs.math.robust_xcorrcoef(yfilt,maxlags,mask);
else
    [R,p]=nirs.math.xcorrcoef(yfilt,maxlags,mask);
end




switch lower(avg_method)
    case 'max'
        dfe = mean(sum(mask)) - 2;

        for i=1:size(R,1)
            for j=1:size(R,2)
                ind = find(abs(R(i,j,:))==max(abs(R(i,j,:))),1);
                R(i,j,1) = R(i,j,ind);
                p(i,j,1) = p(i,j,ind);
            end
        end
        R = R(:,:,1);
        p = p(:,:,1);
    case 'mean'
        dfe = mean(sum(mask))*size(R,3) - 2;

        R = tanh(mean(atanh(R),3));
        R=sign(real(R)).*abs(R);  % the 1's on the diag cause some issues
        Tstat = R .* sqrt((dfe-2) ./ (1 - R.^2));
        p = 2*nirs.math.tpvalue(-abs(Tstat),dfe-2);
        p=p+eye(size(p));
   case 'median'
        dfe = mean(sum(mask))*size(R,3) - 2;

        R = tanh(median(atanh(R),3));
        R=sign(real(R)).*abs(R);  % the 1's on the diag cause some issues
        Tstat = R .* sqrt((dfe-2) ./ (1 - R.^2));
        p = 2*nirs.math.tpvalue(-abs(Tstat),dfe-2);
        p=p+eye(size(p));     
end


