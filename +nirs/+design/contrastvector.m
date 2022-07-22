function [C,haserror] = contrastvector(str,conditions,basis)

%% Examples:
% str = 'stim[3:10]'
% str='stim_channel1[8:20]-stim_channel1[1:3]+stim_channel1[1:2]'

haserror=false;

if(nargin<3)
    basis=[];
end






if(iscellstr(str) | iscell(str))
    for idx=1:length(str)
        [C(idx,:),haserror(idx)] = nirs.design.contrastvector(str{idx},conditions,basis);
    end
    haserror=any(haserror);
    return;
end





% remove spaces
for i=1:length(conditions)
    conditions{i}(strfind(conditions{i},' '))=[];
end
str(strfind(str,' '))=[];

% take care of the simple case
if(ismember(str,conditions))
    C=1*ismember(conditions,str)';
    return;
end

%Parse the str
lst=sort([1 strfind(str,'-') strfind(str,'+') length(str)+1]);
for i=2:length(lst)
    S{i-1}=str(lst(i-1):lst(i)-1);
end

cnt=1;
for i=1:length(S)
    slocal=S{i};
    if(~isempty(strfind(slocal,'x')))
        multiplier{cnt}=str2num(slocal(1:strfind(slocal,'x')-1));
        if(isempty(multiplier{cnt}))
            multiplier{cnt}=1;
        end
    elseif(isempty([strfind(slocal,'-') strfind(slocal,'+')]))
        multiplier{cnt}=1;
    elseif(~isempty(strfind(slocal,'+')))
        multiplier{cnt}=1;
    else
        multiplier{cnt}=-1;
    end
    cond{cnt}=slocal(max([0 strfind(slocal,'-') strfind(slocal,'+') ...
        strfind(slocal,'x')])+1:end);
    cond{cnt}=[cond{cnt}(1:min([strfind(cond{cnt},'[')-1 length(cond{cnt})]))...
        cond{cnt}(min([strfind(cond{cnt},']')+1 length(cond{cnt})+1]):end)];
    
%     
%     for i=1:length(cond)
%         if~(basis.base.iskey(cond{i}))
%             cond{i}=cond{i}(1:min(strfind(cond{i},':'))-1);
%         end
%     end
    
    
    
    if(isempty(strfind(slocal,'[')))
        indices{cnt}=[];
    else
        if(~isempty(strfind(slocal,'canonical')) || (~isempty(strfind(slocal,'gamma'))))
            m=multiplier{cnt};
            if(basis.base.iskey(cond{cnt}))
                base=basis.base(cond{cnt});
            else
                 base=basis.base('default');
            end
            stim=basis.stim(cond{cnt});
            stim.onset=0;
            stim.amp=1;
            dur=max(stim.dur);
            lenHRF=90;
            t=[0:1/basis.Fs:(max(dur)+lenHRF)];
            
            if(~isempty(strfind(slocal,'canonical')))
            [X, names] = nirs.design.createDesignMatrix( Dictionary({cond{cnt}}, {stim}), t, ...
                Dictionary({'default'}, {nirs.design.basis.Canonical}));
            else
                [X, names] = nirs.design.createDesignMatrix( Dictionary({cond{cnt}}, {stim}), t, ...
                Dictionary({'default'}, {nirs.design.basis.Gamma}));
            end
            
            [X2, ~] = nirs.design.createDesignMatrix( Dictionary({cond{cnt}}, {stim}), t, ...
                Dictionary({'default'}, {base}));
            
            samples=[]; cnt2=1; weights=[];
            for i=1:length(conditions)
                if(~isempty(strfind(conditions{i},[cond{cnt} ':'])))
                    samples(cnt2)=str2num(conditions{i}(strfind(conditions{i},':')+1:end));
                    weights(cnt2)=sum(X2(:,samples(cnt2)).*X);
                    cnt2=cnt2+1;
                end
            end
            c=cond{cnt};
            weights=m*weights/max(weights);            
            for i=1:length(weights)
                multiplier{cnt}=weights(i);
                cond{cnt}=c;
                indices{cnt}=samples(i);
                cnt=cnt+1;
            end
            cnt=cnt-1;
        elseif(~isempty(strfind(slocal(strfind(slocal,'['):strfind(slocal,']')),'s')))
            win=slocal(strfind(slocal,'['):strfind(slocal,']'));
            win(strfind(win,'s'))=[];
            win=[win(1:strfind(win,':')) num2str(1/basis.Fs) win(strfind(win,':'):end)];
            win=str2num(win);
            m=multiplier{cnt};
            if(basis.base.iskey(cond{cnt}))
                base=basis.base(cond{cnt});
            else
                 base=basis.base('default');
            end
            stim=basis.stim(cond{cnt});
            stim.onset=0;
            stim.amp=1;
            dur=max(stim.dur);
            lenHRF=90;
            t=[0:1/basis.Fs:(max(dur)+lenHRF)];
            [X2, ~] = nirs.design.createDesignMatrix( Dictionary({cond{cnt}}, {stim}), t, ...
                Dictionary({'default'}, {base}));
            
            X2(t>=max(win),:)=0;
            X2(t<min(win),:)=0;
            weights=sum(X2,1);
           
            weights=m*weights/max(weights);   
             c=cond{cnt};
            for i=1:length(weights)
                multiplier{cnt}=weights(i);
                cond{cnt}=c;
                indices{cnt}=i;
                cnt=cnt+1;
            end
            cnt=cnt-1;
            
         
            
        else
            indices{cnt}=str2num(slocal(strfind(slocal,'['):strfind(slocal,']')));
        end
    end
    cnt=cnt+1;
end


% Remove the numbering off of the stim condition names
for i=1:length(conditions)
    l=strfind(conditions{i},':');
    if(isempty(l))
        l=length(conditions{i});
    else
        l=l-1;
    end
    if(double(conditions{i}(l))>=48 & ...
            double(conditions{i}(l))<=57) %  between 0-9
            names{i}=conditions{i}(1:l-1);
    else
        names{i}=conditions{i};
    end
    block_ind = strfind(conditions{i},'◄');
    if ~isempty(block_ind)
        conditions{i} = conditions{i}(1:block_ind-1);
    end
end
names=unique(names);


C = zeros(1,length(conditions));
% cond{1}(strfind(cond{1},' '))=[];

for idx=1:length(cond)
    if(isempty(indices{idx}))
        lst=find(ismember(conditions,cond{idx}));
        if(length(lst)==0)
            warning(['unable to find name: ' cond{idx}]);
            haserror=true;
        end
    else
        lst=[];
        for i=1:length(indices{idx})
            l=[]; %strfind(cond{idx},':');
            if(isempty(l))
                l=length(cond{idx});
            else
                l=l-1;
            end
            
            s=['00' num2str(indices{idx}(i))];
            lst2=find(ismember(conditions,[cond{idx}(1:l)...
                ':' s(end-1:end) cond{idx}(l+1:end)]));
            if(isempty(lst2))
                warning(['unable to find name: ' cond{idx}(1:l)...
                ':' s(end-1:end) cond{idx}(l+1:end)]);
            haserror=true;
            end
            lst=[lst lst2];
        end
        
    end
   
    C(1,lst)=C(1,lst)+ones(1,length(lst))*multiplier{idx};
end

return

% %% Examples:
% % str = 'stim[3:10]'
% % str='stim_channel1[8:20]-stim_channel1[1:3]+stim_channel1[1:2]'
% 
% if(iscellstr(str) | iscell(str))
%     for idx=1:length(str)
%         C(idx,:) = nirs.design.contrastvector(str{idx},conditions);
%     end
%     return;
% end
% 
% % remove spaces
% for i=1:length(conditions)
%     conditions{i}(strfind(conditions{i},' '))=[];
% end
% str(strfind(str,' '))=[];
% 
% % Sort the variables so that the input order doesn't matter (B:A => A:B) and (D:B-C:B => B:D-B:C)
% sub_parts = strsplit(str,'-');
% for j = 1:length(sub_parts)
%     add_parts = strsplit(sub_parts{j},'+');
%     for k = 1:length(add_parts)
%         terms = strsplit(add_parts{k},':');
%         for l = length(terms):-1:1
%             if all(isstrprop(terms{l},'digit'))
%                 terms{l-1} = [terms{l-1} ':' terms{l}];
%                 terms(l)=[];
%             end
%         end
%         add_parts{k} = strjoin(sort(terms),':');
%     end
%     sub_parts{j} = strjoin(sort(add_parts),'+');
% end
% str = strjoin(sub_parts,'-');
% 
% for i = 1:length(conditions)
%     sub_parts = strsplit(conditions{i},'-');
%     for j = 1:length(sub_parts)
%         add_parts = strsplit(sub_parts{j},'+');
%         for k = 1:length(add_parts)
%             terms = strsplit(add_parts{k},':');
%             for l = length(terms):-1:1
%                 if all(isstrprop(terms{l},'digit'))
%                     terms{l-1} = [terms{l-1} ':' terms{l}];
%                     terms(l)=[];
%                 end
%             end
%             add_parts{k} = strjoin(sort(terms),':');
%         end
%         sub_parts{j} = strjoin(sort(add_parts),'+');
%     end
%     conditions{i} = strjoin(sub_parts,'-');
% end
% 
% % take care of the simple case
% if(ismember(str,conditions))
%     C=1*ismember(conditions,str)';
%     return;
% end
% 
% %Parse the str
% lst=sort([1 strfind(str,'-') strfind(str,'+') length(str)+1]);
% for i=2:length(lst)
%     S{i-1}=str(lst(i-1):lst(i)-1);
% end
% 
% for i=1:length(S)
%     slocal=S{i};
%     if(isempty([strfind(slocal,'-') strfind(slocal,'+')]))
%         multiplier{i}=1;
%     elseif(~isempty(strfind(slocal,'*')))
%         multiplier{i}=str2num(slocal(1:strfind(slocal,'*')-1));
%     elseif(~isempty(strfind(slocal,'+')))
%         multiplier{i}=1;
%     else
%         multiplier{i}=-1;
%     end
%     cond{i}=slocal(max([0 strfind(slocal,'-') strfind(slocal,'+') ...
%         strfind(slocal,'*')])+1:end);
%     cond{i}=[cond{i}(1:min([strfind(cond{i},'[')-1 length(cond{i})]))...
%         cond{i}(min([strfind(cond{i},']')+1 length(cond{i})+1]):end)];
%     if(isempty(strfind(slocal,'[')))
%         indices{i}=[];
%     else
%         indices{i}=str2num(slocal(strfind(slocal,'['):strfind(slocal,']')));
%     end
% end
% 
% 
% % Remove the numbering off of the stim condition names
% for i=1:length(conditions)
%     l=strfind(conditions{i},':');
%     if(isempty(l))
%         l=length(conditions{i});
%     else
%         l=l-1;
%     end
%     if(double(conditions{i}(l))>=48 & ...
%             double(conditions{i}(l))<=57) %  between 0-9
%             names{i}=conditions{i}(1:l-1);
%     else
%         names{i}=conditions{i};
%     end
% end
% names=unique(names);
% 
% 
% C = zeros(1,length(conditions));
% % cond{1}(strfind(cond{1},' '))=[];
% 
% for idx=1:length(cond)
%     if(isempty(indices{idx}))
%         lst=find(ismember(conditions,cond{idx}));
%         if(length(lst)==0)
%             warning(['unable to find name: ' cond{idx}]);
%         end
%     else
%         lst=[];
%         for i=1:length(indices{idx})
%             l=strfind(conditions{idx},':');
%             if(isempty(l))
%                 l=length(conditions{idx});
%             else
%                 l=l-1;
%             end
%             
%             s=['00' num2str(indices{idx}(i))];
%             lst2=find(ismember(conditions,[cond{idx}(1:l)...
%                 ':' s(end-1:end) cond{idx}(l+1:end)]));
%             if(isempty(lst2))
%                 warning(['unable to find name: ' cond{idx}(1:l)...
%                 ':' s(end-1:end) cond{idx}(l+1:end)]);
%             end
%             lst=[lst lst2];
%         end
%         
%     end
%    
%     C(1,lst)=C(1,lst)+ones(1,length(lst))*multiplier{idx};
% end
% 
% return