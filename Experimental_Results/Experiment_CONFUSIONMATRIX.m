%% Experiment Classification of non-convex penalty
clear
DP = [50,300,1350];
density= {'d10','d20'};
MFoldername = strcat('D:\TEMPORAL_THESIS\RESULTS_EXPERIMENT_BIC_',density,'_refit');

n=15;
for p=1:3 % vary lag to see if we choose model lag wrong, what will happen to the solution
    [P,~] = offdiagJSS(n,p,1);
    idx =efficient_vect([n,p,1]);
    indPx = efficient_Px(P,n,p,1);
    for dd=1:2
        MFoldername = strcat('D:\TEMPORAL_THESIS\RESULTS_EXPERIMENT_BIC_',density{dd},'_refit');
        ALL_MAT_ncvx = zeros(10,10);
        ALL_MAT_cvx = zeros(10,10);
        for b=1:10
            for f=1:2
                for dpx = 1:3
                    tic;
                    RR1 = zeros(10,10);
                    RR2 = zeros(10,10);
                    Rname = strcat('DATA_BANK_',int2str(b),'_D',int2str(f),'_',density{dd}); % PASS
                    load(strcat('D:\TEMPORAL_THESIS\DATA_CONFERENCE_CLASSIFICATION\',Rname)) % PASS
                    for mb=1:10
                        Mname_ncvx= strcat(MFoldername,'\ESTIMATED_BANK_',int2str(mb),'_D',int2str(f),'_',int2str(DP(dpx)));
                        Mname_cvx= strcat(MFoldername,'\ESTIMATED_BANK_',int2str(mb),'_D',int2str(f),'_',int2str(DP(dpx)),'_CVX');
                        load(Mname_ncvx)
                        Oncvx(mb)=E(1);
                        load(Mname_cvx)
                        Ocvx(mb) = E(1);
                    end
                    for t=1:20
                        [argmax_bncvx,argmax_bcvx] = CFmat(data.y(:,1:50,1,t),Oncvx,Ocvx,p,P,idx,indPx);
                        RR1(b,argmax_bncvx) = RR1(b,argmax_bncvx)+1;
                        ALL_MAT_ncvx(b,argmax_bncvx)=ALL_MAT_ncvx(b,argmax_bncvx)+1;
                        %                 argmax_bcvx = CFmat(data.y(:,:,1,t),O1,2,P,idx,indPx);
                        RR2(b,argmax_bcvx) = RR2(b,argmax_bcvx)+1;
                        ALL_MAT_cvx(b,argmax_bcvx)=ALL_MAT_cvx(b,argmax_bcvx)+1;
                    end
                    F.ncvx(b,f,dpx).confusion_matrix = RR1;
                    F.cvx(b,f,dpx).confusion_matrix = RR2;
                    
                    toc;
                    
                end
            end
        end
        F.cvxmat = ALL_MAT_cvx;
        F.ncvxmat = ALL_MAT_ncvx;
        
        save(strcat('D:\TEMPORAL_THESIS\CONFUSION_MATRIX_FIX\confusion_matrix_',density{dd},'_p',int2str(p),'_dp50'),'F')
        clear F
    end
    
end
%