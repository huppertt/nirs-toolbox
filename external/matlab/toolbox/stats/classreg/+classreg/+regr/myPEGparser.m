function treestruct = myPEGparser(str,rules)
%*************************************************************
            % Meaning of s and t:
            %
            % treestruct = classreg.regr.LinearMixedFormula.p.parse(f.str); 
            %
            %  s = treestruct.string = Formula string entered by the user
            %  t = treestruct.tree   = parse tree of recorded rules
            %
            %  t is a 5-by-N array of integers. Each column represents a 
            %  node in the tree with rows containing the following info:
            %      rule number for node
            %      start index of match for node
            %      end index of match for node
            %      offset to parent (subtract) of node
            %      number of nodes in the subtree rooted at the node
            %  Any child nodes start directly to the right of the parent
            %  node and appear in order left to right. Thus the entire
            %  subtree rooted at node N is between columns N and
            %  N+TREE(5,N)-1, inclusive.
            %**************************************************************
            
treestruct.str=str;
okchars=[43 45 65:90 97:122 48:57];
okchars2=[65:90 97:122 48:57];

Idx=[1:length(str)];
Idx(strfind(str,' '))=[];
str(strfind(str,' '))=[];


%formula
t(1,1) =rules.Formula;
t(2,1) =Idx(1);
t(3,1) =Idx(end);
t(4,1) =0;

% now the response model
resp=str(1:strfind(str,'~')-1);
idxresp=Idx(1:strfind(str,'~')-1);

t(1,2)=rules.Response;
t(2,2)=idxresp(1);
t(3,2)=idxresp(end);
t(4,2)=1;

lst=[0 find(~ismember(double(resp),okchars)) length(resp)+1];
cnt=3;

nresp=0;
for i=1:length(lst)-1
    t(1,cnt)=rules.ResponseVar;
    t(2,cnt)=idxresp(lst(i)+1);
    t(3,cnt)=idxresp(lst(i+1)-1);
    t(4,cnt)=1;
    t(5,cnt)=1;
    cnt=cnt+1;
    nresp=nresp+1;
end
t(5,2)=nresp+1;

%now the LinearPredictor
linpred=str(strfind(str,'~')+1:end);
linpredIdx=Idx(strfind(str,'~')+1:end);
if(~isempty(strfind(linpred,'(')))
    linpredIdx=linpredIdx(1:strfind(linpred,'(')-1);
    linpred=linpred(1:strfind(linpred,'(')-1);
end

if(ismember(double(linpred(end)),double('-+')))
    linpred(end)=[];
    linpredIdx(end)=[];
end


t(1,cnt)=rules.LinearPredictor;
t(2,cnt)=linpredIdx(1);
t(3,cnt)=linpredIdx(end);
t(4,cnt)=3;
cnt=cnt+1;

t(1,cnt)=rules.Sum;
t(2,cnt)=linpredIdx(1);
t(3,cnt)=linpredIdx(end);
t(4,cnt)=1;
cnt=cnt+1;

npred=0;
lst=[1 find(ismember(double(linpred),double('-+'))) length(linpred)+1];
for i=1:length(lst)-1
    sI=lst(i);
    eI=lst(i+1)-1;
    if(sI>eI)
        lst(i+1)=lst(i);
    else
        npred=npred+1;
        localpred=linpred(sI:eI);
        if(strcmp(linpred(sI),'-'))
            t(1,cnt)=rules.Subend;
        else
            t(1,cnt)=rules.Addend;
        end
        t(2,cnt)=linpredIdx(sI);
        t(3,cnt)=linpredIdx(eI);
        t(4,cnt)=1;
        
        cnt2=cnt;
        npred2=0;
        %TODO sort localpred
        cnt=cnt+1;
        
        if(~isempty(strfind(localpred,':')))
            t(1,cnt)=rules.Interaction;
            t(2,cnt)=linpredIdx(sI);
            t(3,cnt)=linpredIdx(eI);
            t(4,cnt)=1;
            cnt=cnt+1;
             npred=npred+1;
                npred2=npred2+1;
        end
        if(~isempty(strfind(localpred,'*')))
            t(1,cnt)=rules.Product;
            t(2,cnt)=linpredIdx(sI);
            t(3,cnt)=linpredIdx(eI);
            t(4,cnt)=1;
            cnt=cnt+1;
             npred=npred+1;
                npred2=npred2+1;
        end
        if(~isempty(strfind(localpred,'^')))
            t(1,cnt)=rules.Power;
            t(2,cnt)=linpredIdx(sI);
            t(3,cnt)=linpredIdx(eI);
            t(4,cnt)=1;
            cnt=cnt+1;
             npred=npred+1;
                npred2=npred2+1;
        end
        
        lst2=[1 find(ismember(double(localpred),double('-+:*^'))) length(localpred)+1];
        localpredIdx=linpredIdx(sI:eI);
        for j=1:length(lst2)-1
            sI=lst2(j);
            eI=lst2(j+1)-1;
            if(sI>eI)
                lst2(j+1)=lst2(j);
            else
                npred=npred+1;
                npred2=npred2+1;
                localpred2=localpred(sI:eI);
                localpredIdx2=linpredIdx(sI:eI);
                
                localpredIdx2(ismember(double(localpred2),double('-+:*^')))=[];
                localpred2(ismember(double(localpred2),double('-+:*^')))=[];
                
                if(strcmp(localpred(sI),'1'))
                    t(1,cnt)=rules.Intercept;
                
                else
                    t(1,cnt)=rules.PredictorVar;
                end
                t(2,cnt)=localpredIdx2(1);
                t(3,cnt)=localpredIdx2(end);
                t(4,cnt)=1;
                cnt=cnt+1;
            end
        end
        t(5,cnt2)=npred2+1;
    end
end
t(5,find(t(1,:)==rules.Sum))=npred+1;
t(5,find(t(1,:)==rules.LinearPredictor))=npred+2;

%finally, the random effects terms
lst=strfind(str,'(');
for i=1:length(lst)
    RE=str(lst(i):end);
    REidx=Idx(lst(i):end);
    REidx=REidx(1:min(strfind(RE,')')));
    RE=RE(1:min(strfind(RE,')')));
    
    t(1,cnt)=rules.LinearRandomPredictor;
    t(2,cnt)=REidx(1);
    t(3,cnt)=REidx(end);
    t(4,cnt)=11;
    cnt=cnt+1;
    
    eI=strfind(RE,'|');
    t(1,cnt)=rules.LinearPredictor;
    t(2,cnt)=REidx(2);
    t(3,cnt)=REidx(eI-1);
    t(4,cnt)=1;
    cnt=cnt+1;
    
    if(strcmp(RE(2:eI-1),'1'))
        t(1,cnt)=rules.Intercept;
        t(2,cnt)=REidx(2);
        t(3,cnt)=REidx(eI-1);
        t(4,cnt)=1;
        t(5,cnt)=1;
        t(5,cnt-1)=2;
        cnt=cnt+1;
    else
        t(5,cnt-1)=1;
    end
    
    
    t(1,cnt)=rules.GroupingVar;
    t(2,cnt)=REidx(eI+1);
    t(3,cnt)=REidx(end-1);
    t(4,cnt)=3;
    t(5,cnt)=2;
    cnt=cnt+1;
    
    t(1,cnt)=rules.PredictorVar;
    t(2,cnt)=REidx(eI+1);
    t(3,cnt)=REidx(end-1);
    t(4,cnt)=1;
    t(5,cnt)=1;
    cnt=cnt+1;
end

if(isfield(rules,'LinearRandomPredictor'))
    t(5,find(t(1,:)==rules.LinearRandomPredictor))=1+4*length(lst);
end
t(5,1)=size(t,2);


treestruct.tree=t;


