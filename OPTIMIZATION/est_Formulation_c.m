function M = est_Formulation_c(y,A_true,P)
%% This program estimates Ground truth with their original model
% Input are
%       y : 3D array [n,Num,K] which is dimension of timeseries,
%       timepoints, # models respectively
%       A : true VAR parameters [n,n,p,K]
% Output is
%       E : structure containing
%           E.stat
%           E.x_true
%           E.x_est

[n,Num,K] = size(y);
p=2;
gridSize = 60;
Lambda = logspace(-4,0,gridSize);
tmp(K)= struct();
disp('Generating H matrix')
for kk=1:K
    [tmp(kk).H,tmp(kk).Y] = H_gen(y(:,:,kk),p);
end
disp('vectorizing model')
[yc,gc,~] = veccoefmatgroup(reshape([tmp.Y],[n,Num-p,K]),reshape([tmp.H],[n*p,Num-p,K]),A_true,[n,p,K,Num-p]);
disp('calculating Lambda max')
Lmax = lambdamax_grouplasso(gc,yc,[n ,p ,K]);
Lambda = Lambda*Lmax;
Indnorm = max(eig(gc'*gc));
xLS = gc\yc;
indPx = efficient_Px(P,n,p,K);
M.stat.Lambda_max = Lmax;
M.stat.consistency = zeros(2,K,gridSize);
M.stat.bias = zeros(gridSize,1);
for ii=1:gridSize
    [x_APG,~]= nmAPG(gc,yc,indPx,Lambda(ii),2,0.5,[n,p,K],Indnorm,p*K,0.9,0.7,xLS);
    A_est = devect(full(x_APG),n,p,K);
    M.A(:,:,:,:,ii) = A_est;
    [TP,FN,TN,FP] = split_GC_sens(A_true,A_est);
    M.stat.bias(ii) = norm(A_est(:)-A_true(:))/norm(A_true(:));
    M.stat.consistency(:,:,ii) = [TP FN TN FP];
end
end