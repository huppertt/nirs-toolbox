function data=shuffle_demographics(data)
    lstB=1:length(data);
    if(length(data)==1)
        return
    end
    while(true)
        lst=randperm(length(data));
        if(~any(lst-lstB==0))
            break
        end
    end
    data2=data;
    for i=1:length(lst)
        data(lst(i)).demographics=data2(i).demographics;
    end
    
end