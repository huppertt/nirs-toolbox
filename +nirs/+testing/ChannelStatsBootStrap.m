classdef ChannelStatsBootStrap
    
    properties
        formula
        kfold=2;
        grouplevel = true;
        sortfield;
    end
    
    properties (SetAccess = protected)
        Variables
    end
    
    %
    % Different types of ICC [1]
    % Shrout and Fleiss convention	Name in SPSS
    % ICC(1,1)	One-way random single measures
    % ICC(1,k)	One-way random average measures
    % ICC(2,1)	Two-way random single measures (Consistency/Absolute agreement)
    % ICC(2,k)	Two-way random average measures (Consistency/Absolute agreement)
    % ICC(3,1)	Two-way mixed single measures (Consistency/Absolute agreement)
    % ICC(3,k)	Two-way mixed average measures (Consistency/Absolute agreement)
    %
    %
    methods
        function obj = ChannelStatsBootStrap
            obj.formula = 'beta ~ -1 + cond';
            obj.Variables=[];
            obj.sortfield='subject';
        end
        
        
        function obj = run(obj,S)
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            
            subj=unique(demo.(obj.sortfield));
            if(~obj.grouplevel & length(subj)>1)
                % run analysis on a per subject basis
                v=[];
                v2=obj.Variables;
                
                for i=1:length(subj)
                    obj.Variables=[];
                    lst=find(ismember(demo.(obj.sortfield),subj{i}));
                    disp(['processing subject ' subj{i}]);
                    obj=obj.run(S(lst));
                    v=[v; obj.Variables];
                end
                obj.Variables=[v2 v];
                return
            end
            
            % center numeric variables
            n = demo.Properties.VariableNames;
            for i = 1:length(n)
                if all( isnumeric( demo.(n{i}) ) )
                    demo.(n{i}) = demo.(n{i}) - nanmean( demo.(n{i}) );
                end
            end
            
            % preallocate group stats
            G = nirs.core.ChannelStats();
            
            %% loop through files
            W = sparse([]);
            iW = sparse([]);
            
            b = [];
            vars = table();
            for i = 1:length(S)
                % coefs
                if ~isempty(strfind(obj.formula(1:strfind(obj.formula,'~')-1),'tstat'))
                    b = [b; S(i).tstat];
                else
                    b = [b; S(i).beta];
                end
                
                % whitening transform
                
                
                
                [u, s, ~] = svd(S(i).covb, 'econ');
                %W = blkdiag(W, diag(1./diag(sqrt(s))) * u');
                W = blkdiag(W, pinv(s).^.5 * u');
                
                iW = blkdiag(iW, u*sqrt(s) );
                
                
                %                L = chol(S(i).covb,'upper');
                %                W = blkdiag(W,pinv(L));
                
                % table of variables
                file_idx = repmat(i, [size(S(i).beta,1) 1]);
                
                if(~isempty(demo))
                    vars = [vars;
                        [table(file_idx) S(i).variables repmat(demo(i,:), [size(S(i).beta,1) 1])]
                        ];
                else
                    vars = [vars; ...
                        [table(file_idx) S(i).variables]];
                end
            end
            
            % sort
            [vars, idx] = nirs.util.sortrows(vars, {'source', 'detector', 'type'});
            
            % list for first source
            [sd, ~,lst] = nirs.util.uniquerows(table(vars.source, vars.detector, vars.type));
            sd.Properties.VariableNames = {'source', 'detector', 'type'};
            
            %% design mats
            for c = 1:height(vars)
                block_ind = strfind(vars.cond{c},'◄');
                if ~isempty(block_ind)
                    vars.cond{c} = vars.cond{c}(1:block_ind-2);
                end
            end
            
            tmp = vars(lst == 1, :);
            
            beta = randn(size(tmp,1), 1);
            
            nRE=max(1,length(strfind(obj.formula,'|')));
            warning('off','stats:LinearMixedModel:IgnoreCovariancePattern');
            
            obj.formula=nirs.util.verify_formula([table(beta) tmp], obj.formula,true);
            respvar = obj.formula(1:strfind(obj.formula,'~')-1);
            
            lm1 = fitlme([table(beta,'VariableNames',{respvar}) tmp], obj.formula, 'dummyVarCoding',...
                'full', 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));
            
            tb2=[table(beta,'VariableNames',{respvar}) tmp];
            
            tb=lm1.VariableInfo;
            names=tb(tb.InModel & tb.IsCategorical,:).Properties.RowNames;
            names={names{:},'type'};
            v=tb2(:,ismember(tb2.Properties.VariableNames,names));
            
            a=v{:,1};
            for i=2:size(v,2)
                a=strcat(a,v{:,i});
            end
            cv=cvpartition(a,'kfold',obj.kfold);
            X = lm1.designMatrix('Fixed');
            Z = lm1.designMatrix('Random');
            
            nchan = max(lst);
            
            X = kron(speye(nchan), X);
            Z = kron(speye(nchan), Z);
            
            
            %% put them back in the original order
            vars(idx,:) = vars;
            X(idx, :)   = X;
            Z(idx, :)   = Z;
            beta        = b; % already in correct order
            
            for i=1:obj.kfold
                lstTr=find(kron(ones(nchan,1), cv.test(i)));
                Coef(:,i) = nirs.math.fitlme(X(lstTr,:),beta(lstTr),Z(lstTr,:),false,false,false);
                
            end
            
            obj.Variables=[obj.Variables Coef];
            
        end
        
        function disp(obj)
            disp(['NIRS Bootstrap model (' obj.formula ')']);
            if(~obj.grouplevel)
                disp('subject level results');
            end
            disp([num2str(obj.kfold) '-fold testing']);
            disp('-----------------------');
            [r, LB, UB,~,~,~, p] = ICC(obj,'1-1');
            disp(['ICC(1,1) One-way random single ']);
            disp(['                               absolute ' num2str(r) ' [' num2str(LB) '-' num2str(UB) '] p=' num2str(p)]);
               [r, LB, UB,~,~,~, p] = ICC(obj,'C-1');
            disp(['                            consistency ' num2str(r) ' [' num2str(LB) '-' num2str(UB) '] p=' num2str(p)]);          
            [r, LB, UB,~,~,~, p] = ICC(obj,'A-1');
            disp(['                            reliability ' num2str(r) ' [' num2str(LB) '-' num2str(UB) '] p=' num2str(p)]);  
            
            
            [r, LB, UB,~,~,~, p] = ICC(obj,'1-k');
            disp(['ICC(1,k) One-way random average ']);
            disp(['                              absolute ' num2str(r) ' [' num2str(LB) '-' num2str(UB) '] p=' num2str(p)]);
            [r, LB, UB,~,~,~, p] = ICC(obj,'C-k');
                 [r, LB, UB,~,~,~, p] = ICC(obj,'C-1');
            disp(['                            consistency ' num2str(r) ' [' num2str(LB) '-' num2str(UB) '] p=' num2str(p)]);          
            [r, LB, UB,~,~,~, p] = ICC(obj,'A-k');
            disp(['                            reliability ' num2str(r) ' [' num2str(LB) '-' num2str(UB) '] p=' num2str(p)]);
            

            
            
        end
        
        function [r, LB, UB, F, df1, df2, p] = ICC(obj, type, alpha, r0)
            M=obj.Variables;
            % Intraclass correlation
            %   [r, LB, UB, F, df1, df2, p] = ICC(M, type, alpha, r0)
            %
            %   M is matrix of observations. Each row is an object of measurement and
            %   each column is a judge or measurement.
            %
            %   'type' is a string that can be one of the six possible codes for the desired
            %   type of ICC:
            %       '1-1': The degree of absolute agreement among measurements made on
            %         randomly seleted objects. It estimates the correlation of any two
            %         measurements.
            %       '1-k': The degree of absolute agreement of measurements that are
            %         averages of k independent measurements on randomly selected
            %         objects.
            %       'C-1': case 2: The degree of consistency among measurements. Also known
            %         as norm-referenced reliability and as Winer's adjustment for
            %         anchor points. case 3: The degree of consistency among measurements maded under
            %         the fixed levels of the column factor. This ICC estimates the
            %         corrlation of any two measurements, but when interaction is
            %         present, it underestimates reliability.
            %       'C-k': case 2: The degree of consistency for measurements that are
            %         averages of k independent measurements on randomly selected
            %         onbjectgs. Known as Cronbach's alpha in psychometrics. case 3:
            %         The degree of consistency for averages of k independent
            %         measures made under the fixed levels of column factor.
            %       'A-1': case 2: The degree of absolute agreement among measurements. Also
            %         known as criterion-referenced reliability. case 3: The absolute
            %         agreement of measurements made under the fixed levels of the column factor.
            %       'A-k': case 2: The degree of absolute agreement for measurements that are
            %         averages of k independent measurements on randomly selected objects.
            %         case 3: he degree of absolute agreement for measurements that are
            %         based on k independent measurements maded under the fixed levels
            %         of the column factor.
            %
            %       ICC is the estimated intraclass correlation. LB and UB are upper
            %       and lower bounds of the ICC with alpha level of significance.
            %
            %       In addition to estimation of ICC, a hypothesis test is performed
            %       with the null hypothesis that ICC = r0. The F value, degrees of
            %       freedom and the corresponding p-value of the this test are
            %       reported.
            %
            %       (c) Arash Salarian, 2008
            %
            %       Reference: McGraw, K. O., Wong, S. P., "Forming Inferences About
            %       Some Intraclass Correlation Coefficients", Psychological Methods,
            %       Vol. 1, No. 1, pp. 30-46, 1996
            %
            if(nargin<2)
                type='1-1';
            end
            if nargin < 3
                alpha = .05;
            end
            
            if nargin < 4
                r0 = 0;
            end
            
            [n, k] = size(M);
            
            SStotal = var(M(:)) *(n*k - 1);
            MSR = var(mean(M, 2)) * k;
            MSW = sum(var(M,0, 2)) / n;
            MSC = var(mean(M, 1)) * n;
            MSE = (SStotal - MSR *(n - 1) - MSC * (k -1))/ ((n - 1) * (k - 1));
            
            switch type
                case '1-1'
                    [r, LB, UB, F, df1, df2, p] = ICC_case_1_1(MSR, MSE, MSC, MSW, alpha, r0, n, k);
                case '1-k'
                    [r, LB, UB, F, df1, df2, p] = ICC_case_1_k(MSR, MSE, MSC, MSW, alpha, r0, n, k);
                case 'C-1'
                    [r, LB, UB, F, df1, df2, p] = ICC_case_C_1(MSR, MSE, MSC, MSW, alpha, r0, n, k);
                case 'C-k'
                    [r, LB, UB, F, df1, df2, p] = ICC_case_C_k(MSR, MSE, MSC, MSW, alpha, r0, n, k);
                case 'A-1'
                    [r, LB, UB, F, df1, df2, p] = ICC_case_A_1(MSR, MSE, MSC, MSW, alpha, r0, n, k);
                case 'A-k'
                    [r, LB, UB, F, df1, df2, p] = ICC_case_A_k(MSR, MSE, MSC, MSW, alpha, r0, n, k);
            end
            
            
        end
        
    end
end


%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_1_1(MSR, MSE, MSC, MSW, alpha, r0, n, k)
r = (MSR - MSW) / (MSR + (k-1)*MSW);

F = (MSR/MSW) * (1-r0)/(1+(k-1)*r0);
df1 = n-1;
df2 = n*(k-1);
p = 1-fcdf(F, df1, df2);

FL = (MSR/MSW) / finv(1-alpha/2, n-1, n*(k-1));
FU = (MSR/MSW) * finv(1-alpha/2, n*(k-1), n-1);

LB = (FL - 1) / (FL + (k-1));
UB = (FU - 1) / (FU + (k-1));
end

%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_1_k(MSR, MSE, MSC, MSW, alpha, r0, n, k)
r = (MSR - MSW) / MSR;

F = (MSR/MSW) * (1-r0);
df1 = n-1;
df2 = n*(k-1);
p = 1-fcdf(F, df1, df2);

FL = (MSR/MSW) / finv(1-alpha/2, n-1, n*(k-1));
FU = (MSR/MSW) * finv(1-alpha/2, n*(k-1), n-1);

LB = 1 - 1 / FL;
UB = 1 - 1 / FU;
end
%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_C_1(MSR, MSE, MSC, MSW, alpha, r0, n, k)
r = (MSR - MSE) / (MSR + (k-1)*MSE);

F = (MSR/MSE) * (1-r0)/(1+(k-1)*r0);
df1 = n - 1;
df2 = (n-1)*(k-1);
p = 1-fcdf(F, df1, df2);

FL = (MSR/MSE) / finv(1-alpha/2, n-1, (n-1)*(k-1));
FU = (MSR/MSE) * finv(1-alpha/2, (n-1)*(k-1), n-1);

LB = (FL - 1) / (FL + (k-1));
UB = (FU - 1) / (FU + (k-1));
end

%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_C_k(MSR, MSE, MSC, MSW, alpha, r0, n, k)
r = (MSR - MSE) / MSR;

F = (MSR/MSE) * (1-r0);
df1 = n - 1;
df2 = (n-1)*(k-1);
p = 1-fcdf(F, df1, df2);

FL = (MSR/MSE) / finv(1-alpha/2, n-1, (n-1)*(k-1));
FU = (MSR/MSE) * finv(1-alpha/2, (n-1)*(k-1), n-1);

LB = 1 - 1 / FL;
UB = 1 - 1 / FU;
end
%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_A_1(MSR, MSE, MSC, MSW, alpha, r0, n, k)
r = (MSR - MSE) / (MSR + (k-1)*MSE + k*(MSC-MSE)/n);

a = (k*r0) / (n*(1-r0));
b = 1 + (k*r0*(n-1))/(n*(1-r0));
F = MSR / (a*MSC + b*MSE);
%df2 = (a*MSC + b*MSE)^2/((a*MSC)^2/(k-1) + (b*MSE)^2/((n-1)*(k-1)));

a = k*r/(n*(1-r));
b = 1+k*r*(n-1)/(n*(1-r));
v = (a*MSC + b*MSE)^2/((a*MSC)^2/(k-1) + (b*MSE)^2/((n-1)*(k-1)));

df1 = n - 1;
df2 = v;
p = 1-fcdf(F, df1, df2);

Fs = finv(1-alpha/2, n-1, v);
LB = n*(MSR - Fs*MSE)/(Fs*(k*MSC + (k*n - k - n)*MSE) + n*MSR);

Fs = finv(1-alpha/2, v, n-1);
UB = n*(Fs*MSR-MSE)/(k*MSC + (k*n - k - n)*MSE + n*Fs*MSR);
end

%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_A_k(MSR, MSE, MSC, MSW, alpha, r0, n, k)
r = (MSR - MSE) / (MSR + (MSC-MSE)/n);

c = r0/(n*(1-r0));
d = 1 + (r0*(n-1))/(n*(1-r0));
F = MSR / (c*MSC + d*MSE);
%df2 = (c*MSC + d*MSE)^2/((c*MSC)^2/(k-1) + (d*MSE)^2/((n-1)*(k-1)));

a = k*r/(n*(1-r));
b = 1+k*r*(n-1)/(n*(1-r));
v = (a*MSC + b*MSE)^2/((a*MSC)^2/(k-1) + (b*MSE)^2/((n-1)*(k-1)));

df1 = n - 1;
df2 = v;
p = 1-fcdf(F, df1, df2);

Fs = finv(1-alpha/2, n-1, v);
LB = n*(MSR - Fs*MSE)/(Fs*(MSC-MSE) + n*MSR);

Fs = finv(1-alpha/2, v, n-1);
UB = n*(Fs*MSR - MSE)/(MSC - MSE + n*Fs*MSR);
end
