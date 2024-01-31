function tbl1020 = register_headsize(headsize,tbl1020)
% This function will transform the 10-20 table to match the head size
% Headsize is a dictionary of the form:
%   headsize=Dictionary();
%   headsize('lpa-Cz-rpa')=346;   - specifies a specific distance (in mm) from two points
%   headsize('Iz-Cz-nas')=373;
%   headsize('circumference')=523  - or the circumfurance
%

if(nargin<2)
    tbl1020=nirs.util.list_1020pts('?');
end

XYZ=[tbl1020.X tbl1020.Y tbl1020.Z];

keys=headsize.keys;
idx1=find(ismember(lower(tbl1020.Name),'lpa'));
idx2=find(ismember(lower(tbl1020.Name),'rpa'));
idx3=find(ismember(lower(tbl1020.Name),'cz'));
idx4=find(ismember(lower(tbl1020.Name),'nas'));
idx5=find(ismember(lower(tbl1020.Name),'iz'));

for i=1:headsize.count
    val=double(headsize(keys{i}));
    str=lower(keys{i});
    str(strfind(str,' '))=[];
        
    if(~isempty(strfind(str,'circ')) || ...
            ~isempty(strfind(str,'circumference')))
            %a=obj.LR_distance/2;
            %b=obj.AP_distance/2;
           
           a=(XYZ(idx1,:)-XYZ(idx2,:));
           b=(XYZ(idx4,:)-XYZ(idx5,:));
          cost{i}=@(s)abs(circumference(.5*norm(s.*a),.5*norm(s.*b))-val); 
    elseif(strcmp(lower(str),'lpa-cz-rpa') || ...
                strcmp(lower(str),'rpa-cz-lpa'))
            
          %a=obj.LR_distance/2;
          %b=obj.IS_distance;
           
          a=(XYZ(idx1,:)-XYZ(idx2,:));
          b=XYZ(idx3,:)-.5*(XYZ(idx1,:)-XYZ(idx2,:));
          cost{i}=@(s)abs(arcdistance(.5*norm(s.*a),norm(s.*b))-val);
          
    elseif(strcmp(lower(str),'nas-cz-iz') || ...
                strcmp(lower(str),'iz-cz-nas'))
            
          %a=obj.AP_distance/2;
          %b=obj.IS_distance;
           
          a=(XYZ(idx4,:)-XYZ(idx5,:));
          XYZ(idx3,:)-.5*(XYZ(idx1,:)-XYZ(idx2,:));
          cost{i}=@(s)abs(arcdistance(.5*norm(s.*a),norm(s.*b))-val);
    else
        warning('head size key not recognized');
    end        
end


if(~isempty(ver('optim')))
    opt = optimoptions('lsqnonlin', 'MaxFunEvals', 1000,'Display','off');
   
    if(headsize.count>2)
        cstfcn=@(x)vertcat(cell2mat(cellfun(@(a){a(x)},cost)'));
        s=lsqnonlin(cstfcn,[1 1 1],[.5 .5 .5],[2 2 2],opt);
    else
        cstfcn=@(x)vertcat(cell2mat(cellfun(@(a){a(x*[1 1 1])},cost)'));
        s=lsqnonlin(cstfcn,1,.5,2,opt);
        s=s*[1 1 1];
    end
else
    opt = optimset('MaxFunEvals', 1000,'Display','off');
     if(headsize.count>2)
        cstfcn=@(x)max(reshape(vertcat(cell2mat(cellfun(@(a){a(x)},cost)')),[],1));
        %s=lsqnonlin(cstfcn,[1 1 1],[.5 .5 .5],[2 2 2],opt);
        s=fminsearch(cstfcn,[1 1 1],opt);
    else
        cstfcn=@(x)max(reshape(vertcat(cell2mat(cellfun(@(a){a(x*[1 1 1])},cost)')),[],1));
        s=fminbnd(cstfcn,.5,2,opt);
        s=s*[1 1 1];
    end
    
end


tbl1020.X=XYZ(:,1)*s(1);
tbl1020.Y=XYZ(:,2)*s(2);
tbl1020.Z=XYZ(:,3)*s(3);

end

function arclength=arcdistance(a,b)
% This function returns the arc length between two points
  
   arclength=pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)))/2;
end

function circ = circumference(a,b)
   % This finds the head circumference at 10% up
   %a=obj.LR_distance/2;
   %b=obj.AP_distance/2;      
   
    circ=.9*pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)));
end
