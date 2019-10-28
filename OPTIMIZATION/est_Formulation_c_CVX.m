function M = est_Formulation_c_CVX(y,A_true_in,P)
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
A_true = A_true_in(:,:,:,1:K);
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
% Indnorm = max(eig(gc'*gc));
xLS = gc\yc;
% indPx = efficient_Px(P,n,p,K);
M.stat.Lambda_max = Lmax;
M.stat.consistency = zeros(2,4,gridSize);
M.stat.bias = zeros(gridSize,1);
for ii=1:gridSize
    [x,~, ~] = grouplasso_sharedsp(gc, yc, P,0, Lambda(ii), [n,p,K], 10*max([Lambda(ii),10]),Num-p,xLS);
    A_est = devect(full(x),n,p,K);
    M.A(:,:,:,:,ii) = A_est;
    [TP,FN,TN,FP] = split_GC_sens_K(A_true,A_est);
    M.stat.bias(ii) = norm(A_est(:)-A_true(:))/norm(A_true(:));
    M.stat.consistency(:,:,ii) = [TP FN TN FP];
end
end