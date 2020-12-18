function data = permutedata(data,perservespace)

if(nargin<2)
    perservespace=false;
end

for i=1:length(data)
    [yfilt,f] = nirs.math.innovations(data(i).data,fix(data(i).Fs*4));
    
    n=nnz(yfilt<0)/nnz(yfilt);
    if(perservespace)
        yfilt=abs(yfilt).*sign(rand(size(yfilt,1),1)*ones(1,size(yfilt,2))-n);
    else
        yfilt=abs(yfilt).*sign(rand(size(yfilt))-n);
    end
    
       
    for j=1:length(f)
        data(i).data(:,j)=filter(1,f{j},yfilt(:,j));
    end
    
    
end
