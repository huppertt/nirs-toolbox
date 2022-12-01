function Fisher = statecmmvnrfish(Data, Design, Covar, Method, MatrixFormat,CholCovar,isdiagonal)
%STATECMMVNRFISH Fisher information for multivariate normal regression model.
%    Fisher information matrix based on current maximum likelihood or
%    least-squares parameter estimates that account for missing data.

%    Copyright 2006-2007 The MathWorks, Inc.


[NumSamples, NumSeries] = size(Data);
if iscell(Design)
    if (numel(Design) == 1)
        SingleDesign = true;
    else
        SingleDesign = false;
    end
    NumParams = size(Design{1},1);
else
    SingleDesign = false;
    NumParams = size(Design,2);
end

% Step 2 - initialization
if strcmpi(MatrixFormat,'PARAMONLY')
    TotalParams = NumParams;
elseif isdiagonal
    TotalParams = NumParams + NumSeries;
else
    TotalParams = NumParams + (NumSeries * (NumSeries + 1))/2;
end

Fisher = zeros(TotalParams,TotalParams);

if isdiagonal
    % Get diagonal inverse covariance and the indices of the corresponding
    % diagonal elements in the Fisher matrix
    InvCovar = diag(1 ./ diag(Covar));
    diagind = (NumParams*TotalParams+NumParams+1:TotalParams+1:TotalParams^2)';
else
    InvCovar = inv(Covar);
end

% Step 3 - calculate fisher information matrix (not hessian)
if strcmpi(Method,'FISHER')
    
    % Step 4 - do partials wrt Mean
    Count = 0;
    TestMatrix = zeros(NumParams,NumParams);
    if iscell(Design)
        if SingleDesign
            A = CholCovar' \ Design{1};
            for k = 1:NumSamples
                TestMatrix = TestMatrix + A'*A;
            end
            Count = NumSamples;
        else
            for k = 1:NumSamples
                if ~all(isnan(Data(k,:)))
                    Count = Count + 1;
                    A = CholCovar' \ Design{k}';
                    TestMatrix = TestMatrix + A'*A;
                end
            end
        end
    else
        for k = 1:NumSamples
            if ~all(isnan(Data(k,:)))
                Count = Count + 1;
                A = CholCovar' \ Design(k,:);
                TestMatrix = TestMatrix + A'*A;
            end
        end
    end
    TestMatrix = (1.0/Count) .* TestMatrix;
    Fisher(1:NumParams,1:NumParams) = TestMatrix;

    % Step 5 - do partials wrt Covar

    if strcmpi(MatrixFormat,'FULL')
        if isdiagonal
            cdiag = diag(InvCovar);    % inverse covariance diagonal
            fdiag = 0.5 * cdiag.^2;    % Fisher matrix diagonal
            Fisher(diagind) = fdiag;
        else
            GradC1 = zeros(NumSeries,NumSeries);
            GradC2 = zeros(NumSeries,NumSeries);

            i = NumParams;
            for i1 = 1:NumSeries
                for j1 = 1:i1
                    i = i + 1;

                    GradC1(i1,j1) = 1.0;                    % do dC/dtheta(i)
                    GradC1(j1,i1) = 1.0;
                    Temp1 = InvCovar*GradC1;
                    GradC1(i1,j1) = 0.0;                    % undo dC/dtheta(i)
                    GradC1(j1,i1) = 0.0;

                    j = NumParams;
                    for i2 = 1:NumSeries
                        for j2 = 1:i2
                            j = j + 1;

                            if (j <= i)
                                GradC2(i2,j2) = 1.0;        % do dC/dtheta(j)
                                GradC2(j2,i2) = 1.0;

                                Temp2 = InvCovar*GradC2;

                                Fisher(i,j) = 0.5*trace(Temp1*Temp2);
                                Fisher(j,i) = Fisher(i,j);

                                GradC2(i2,j2) = 0.0;        % undo dC/dtheta(j)
                                GradC2(j2,i2) = 0.0;
                            end
                        end
                    end
                end
            end
        end
    end

    Fisher = Count * Fisher;

    % Step 6 - calculate hessian (not fisher information matrix)

else
    % Step 7 - main loop over data records
    Count = 0;

    for kk = 1:NumSamples

        % Step 8 - determine and map available data in current record

        Map = ~isnan(Data(kk,:));
        Available = sum(Map);

        if Available > 0                        % skip over empty records
            Count = Count + 1;

            % Step 9 - construct covariance matrix subarrays

            if iscell(Design)
                if SingleDesign
                    SubDesign = Design{1};
                else
                    SubDesign = Design{kk};
                end
            else
                SubDesign = Design(kk,:);
            end

            if Available < NumSeries
                SubDesign = SubDesign(Map,:);
				if isdiagonal
					InvSubCovar = inv(Covar(Map,Map));
				else
					InvSubCovar = diag(1 ./ diag(Covar(Map,Map)));
				end
            else
                InvSubCovar = InvCovar;
            end

            % Step 10 - do partials wrt Mean for current data record

            TempMatrix = SubDesign' * InvSubCovar * SubDesign;
            Fisher(1:NumParams,1:NumParams) = ...
                Fisher(1:NumParams,1:NumParams) + TempMatrix;
            
            % Step 11 - do partials wrt Covar for current data record

            if strcmpi(MatrixFormat,'FULL')
                if isdiagonal
                    cdiag = diag(InvCovar);       % inverse covariance diag
                    fdiag = 0.5 * cdiag(Map).^2;  % Fisher matrix diagonal
                    Fisher(diagind(Map)) = Fisher(diagind(Map)) + fdiag;
                else
                    GradC1 = zeros(Available,Available);
                    GradC2 = zeros(Available,Available);

                    p1 = 0;
                    i = NumParams;
                    for i1 = 1:NumSeries
                        if Map(i1)
                            p1 = p1 + 1;
                        end

                        q1 = 0;
                        for j1 = 1:i1
                            i = i + 1;

                            if Map(j1)
                                q1 = q1 + 1;
                            end

                            if Map(i1) && Map(j1)
                                GradC1(p1,q1) = 1.0;
                                GradC1(q1,p1) = 1.0;
                                Temp1 = InvSubCovar*GradC1;
                                GradC1(p1,q1) = 0.0;
                                GradC1(q1,p1) = 0.0;
                            end

                            p2 = 0;
                            j = NumParams;
                            for i2 = 1:NumSeries
                                if Map(i2)
                                    p2 = p2 + 1;
                                end

                                q2 = 0;
                                for j2 = 1:i2
                                    j = j + 1;

                                    if Map(j2)
                                        q2 = q2 + 1;
                                    end

                                    % dC/dtheta(i) = dC/dC(i1,j1)
                                    % dC/dtheta(j) = dC/dC(i2,j2)

                                    if (j <= i) && Map(i1) && Map(j1) && Map(i2) && Map(j2)
                                        GradC2(p2,q2) = 1.0;
                                        GradC2(q2,p2) = 1.0;
                                        Temp2 = InvSubCovar*GradC2;
                                        GradC2(p2,q2) = 0.0;
                                        GradC2(q2,p2) = 0.0;

                                        Fisher(i,j) = Fisher(i,j) + 0.5*trace(Temp1*Temp2);
                                        Fisher(j,i) = Fisher(i,j);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
