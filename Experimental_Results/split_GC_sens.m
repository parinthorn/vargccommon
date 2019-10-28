function [TP,FN,TN,FP] = split_GC_sens(A_true,A_est)
% A = [n,n,p,4]; ONLY FOR K=4, Can be extend using for loop
Dim = size(A_true);
GC_true_K = squeeze(A_true(:,:,1,:)) ~=0 ;
GC_est_K = squeeze(A_est(:,:,1,:)) ~=0 ;

ind_common_true = (GC_true_K(:,:,1)~=0) & (GC_true_K(:,:,2)~=0) & ...
                  (GC_true_K(:,:,3)~=0)&  (GC_true_K(:,:,4)~=0);

ind_common_est = (GC_est_K(:,:,1)~=0) & (GC_est_K(:,:,2)~=0) & ...
                  (GC_est_K(:,:,3)~=0)&  (GC_est_K(:,:,4)~=0);

TP_common = sum((ind_common_true) & (ind_common_est),'all')-Dim(1);
TN_common = sum((~ind_common_true) & (~ind_common_est),'all');
FP_common = sum((~ind_common_true) & (ind_common_est),'all');
FN_common = sum((ind_common_true) & (~ind_common_est),'all');
TP_diff = 0;
TN_diff = 0;
FP_diff = 0;
FN_diff = 0;
for kk=1:4
    GC_true = GC_true_K(:,:,kk);
    GC_est = GC_est_K(:,:,kk);
    TP_diff = TP_diff+sum((~ind_common_true & GC_true) & (~ind_common_est & GC_est),'all');
    TN_diff = TN_diff+sum(~(~ind_common_true & GC_true) & ~(~ind_common_est & GC_est),'all');
    FP_diff = FP_diff+sum(~(~ind_common_true & GC_true) & (~ind_common_est & GC_est),'all');
    FN_diff = FN_diff+sum((~ind_common_true & GC_true) & ~(~ind_common_est & GC_est),'all');
end
TP = [TP_common;TP_diff];
TN = [TN_common;TN_diff];
FP = [FP_common;FP_diff];
FN = [FN_common;FN_diff];
end