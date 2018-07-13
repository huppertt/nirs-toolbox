function PSF = inverse_pointspread(Jacob,mesh)
% This function computes the point-spread function from a jacobian
%
% Use:
% mesh.draw(PSF); % to draw

if(~iscell(Jacob)); Jacob={Jacob}; end;

cnt=1;
for j=1:length(Jacob)
    flds=fields(Jacob{j});
    for i=1:length(flds)
        J{cnt}=abs(Jacob{j}.(flds{i}));
        J{cnt}=J{cnt}/normest(J{cnt});
        cnt=cnt+1;
    end
end
J=vertcat(J{:});

% Point spread
[~,~,V]=nirs.math.mysvd(J);

% J = U*S*V'
% inv(J) = V*inv(S)*S' = V*inv(US)
% 
% L=US{1};
% for i=2:length(US)
%     %???? TODO
%     L=[L; US{i}];
% end

%PSF 
% PSF = V* pinv(L)*L*V';  = V*V'  (so we just need to keep V)

% Now compute the distances
disp('computing distance matrix');
X=mesh.nodes(:,1);
X=repmat(X,1,length(X));
X=X-X';
Y=mesh.nodes(:,2);
Y=repmat(Y,1,length(Y));
Y=Y-Y';
Z=mesh.nodes(:,3);
Z=repmat(Z,1,length(Z));
Z=Z-Z';
DIST = sqrt(X.^2+Y.^2+Z.^2);
clear X Y Z;

PSF = abs(V*V');
PSF=PSF./(ones(size(PSF,1),1)*max(PSF,[],1));

PSF=(PSF>1/exp(1));  % FWHM;
PSF = PSF.*DIST;
PSF=sort(PSF,1,'descend');
PSF=median(PSF(1:10,:),1)';



