%% BIC assignment to MODEL in SIMULATION DATASET
% NOTE THAT THE SAVE IS COMMENTED
clear 
clc
Fdataname = 'G:\Shared drives\MASTER_DRIVE\THESIS\THESIS_DATA_conference';
Fmodelnamed10 = 'G:\Shared drives\MASTER_DRIVE\THESIS\RESULTS\FormulationC_Model16Oct19_d10';
Fmodelnamed20 = 'G:\Shared drives\MASTER_DRIVE\THESIS\RESULTS\FormulationC_Model16Oct19_d20';
S1Folder = 'G:\Shared drives\MASTER_DRIVE\THESIS\RESULTS_EXPERIMENT_BIC_d10_refit';
S2Folder = 'G:\Shared drives\MASTER_DRIVE\THESIS\RESULTS_EXPERIMENT_BIC_d20_refit';
n=15;
p=2;
K=4;
[P,~] = offdiagJSS(15,2,4);
indPx = efficient_Px(P,15,2,4);
idx = efficient_vect([15,2,4]);
DP = [50,300,1350];
for b=1:10
    disp(b)
    for f=1:2
        for dpx=1:3
            % LOAD, Data, x_est
            Rname = strcat('\ESTIMATED_BANK_',int2str(b),'_D',int2str(f),'_',int2str(DP(dpx)),'_CVX.mat');
            D10name = strcat('\DATA_BANK_',int2str(b),'_D',int2str(f),'_d10.mat');
            D20name = strcat('\DATA_BANK_',int2str(b),'_D',int2str(f),'_d20.mat');
            load(strcat(Fdataname,D10name))
            yd10 = data.y(:,1:DP(dpx),:,:);
            load(strcat(Fdataname,D20name))
            yd20 = data.y(:,1:DP(dpx),:,:);
            load(strcat(Fmodelnamed10,Rname))
            Ed10 = E;
            load(strcat(Fmodelnamed20,Rname))
            Ed20 = E;
            tic
            for t = 1:4
                A_D10_trial = Ed10(t).M.A;
                A_D20_trial = Ed20(t).M.A;
                dataD10 = yd10(:,:,:,t);
                dataD20 = yd20(:,:,:,t);
                clear tmp1 tmp2
                tmp1(4)= struct();
                tmp2(4)= struct();
                for kk=1:4
                    [tmp1(kk).H,tmp1(kk).Y] = H_gen(dataD10(:,:,kk),2);
                    [tmp2(kk).H,tmp2(kk).Y] = H_gen(dataD20(:,:,kk),2);
                end
                [yc1,gc1,~] = veccoefmatgroup(reshape([tmp1.Y],[n,DP(dpx)-p,K]),reshape([tmp1.H],[n*p,DP(dpx)-p,K]),zeros(n,n,p,K),[n,p,K,DP(dpx)-2]);
                [yc2,gc2,~] = veccoefmatgroup(reshape([tmp2.Y],[n,DP(dpx)-p,K]),reshape([tmp2.H],[n*p,DP(dpx)-p,K]),zeros(n,n,p,K),[n,p,K,DP(dpx)-2]);

                xLS1 = gc1\yc1;
                xLS2 = gc2\yc2;
                for LL=1:60 % Lambda gives different BIC
                    A1tmp = A_D10_trial(:,:,:,:,LL);
                    A2tmp = A_D20_trial(:,:,:,:,LL);
                    x1 = refit(gc1,yc1,P,A1tmp(idx(indPx)),1);
                    x2 = refit(gc2,yc2,P,A2tmp(idx(indPx)),1);
                    Ed10(t).M.A(:,:,:,:,LL) = devect(x1,n,p,K);
                    Ed10(t).M.stat.bic(LL) = BIC(dataD10,Ed10(t).M.A(:,:,:,:,LL));
                    Ed10(t).M.stat.bic_yuan(LL) = BIC_yuan(dataD10,xLS1,Ed10(t).M.A(:,:,:,:,LL),idx);
                    
                    Ed20(t).M.A(:,:,:,:,LL) = devect(x2,n,p,K);
                    Ed20(t).M.stat.bic(LL) = BIC(dataD20,Ed20(t).M.A(:,:,:,:,LL));
                    Ed20(t).M.stat.bic_yuan(LL)= BIC_yuan(dataD20,xLS2,Ed20(t).M.A(:,:,:,:,LL),idx);
                end
                [~,bb] = min(Ed10(t).M.stat.bic);
                Ed10(t).M.stat.argmin_bic = bb;
                [~,bb] = min(Ed20(t).M.stat.bic);
                Ed20(t).M.stat.argmin_bic = bb;                
                [~,bb] = min(Ed10(t).M.stat.bic_yuan);
                Ed10(t).M.stat.argmin_bic_yuan = bb;
                [~,bb] = min(Ed20(t).M.stat.bic_yuan);
                Ed20(t).M.stat.argmin_bic_yuan = bb;    
            end
            toc
%             E = Ed10;
%             save(strcat(S1Folder,Rname),'E')
%             E = Ed20;
%             save(strcat(S2Folder,Rname),'E')
        end
    end
    
    
end

