% This experiment shows that non-convex model is more robust to density changing
% than convex model

clear
clc
S1Folder = 'G:\Shared drives\MASTER_DRIVE\THESIS\RESULTS_EXPERIMENT_BIC_d10';
S2Folder = 'G:\Shared drives\MASTER_DRIVE\THESIS\RESULTS_EXPERIMENT_BIC_d20';
DP= [50 300 1350];
% BIC_PERFORMANCE = zeros(10,2,3,4);
% BIC_YUAN_PERFORMANCE = zeros(10,2,3,4); BIC_yuan
% BIC_PERFORMANCE(10,2,3).bic = struct();
% BIC_PERFORMANCE(10,2,3).bic_yuan = struct();
BIC_avg = zeros(2,4);
BIC_yuan_avg = zeros(2,4);
TYPE = {'','_CVX'};
for c=1:2
for b=1:10
    for f=1:2
        for dpx = 1:3
            Rname = strcat('\ESTIMATED_BANK_',int2str(b),'_D',int2str(f),'_',int2str(DP(dpx)),TYPE{c},'.mat');
            load(strcat(S2Folder,Rname))
            BIC_PERFORMANCE(b,f,dpx,c).bic = zeros(1,4);
            BIC_PERFORMANCE(b,f,dpx,c).bic_yuan = zeros(1,4);
            for t=1:4
                BIC_PERFORMANCE(b,f,dpx,c).bic = BIC_PERFORMANCE(b,f,dpx,c).bic  + E(t).M.stat.consistency(1,:,E(t).M.stat.argmin_bic);
                BIC_PERFORMANCE(b,f,dpx,c).bic_yuan = BIC_PERFORMANCE(b,f,dpx,c).bic_yuan + E(t).M.stat.consistency(1,:,E(t).M.stat.argmin_bic_yuan);
                BIC_avg(c,:) = BIC_avg(c,:)+BIC_PERFORMANCE(b,f,dpx,c).bic;
                BIC_yuan_avg(c,:) = BIC_yuan_avg(c,:)+BIC_PERFORMANCE(b,f,dpx,c).bic_yuan;
            end
        end
    end
end
end

%%

TPRbic = BIC_avg(:,1)./(BIC_avg(:,1)+BIC_avg(:,2));
FPRbic = BIC_avg(:,4)./(BIC_avg(:,3)+BIC_avg(:,4));

TPRbicy = BIC_yuan_avg(:,1)./(BIC_yuan_avg(:,1)+BIC_yuan_avg(:,2));
FPRbicy = BIC_yuan_avg(:,4)./(BIC_yuan_avg(:,3)+BIC_yuan_avg(:,4));

%% Is it invariant under topology ?


for c=1:2
for b=1:10
    for f=1:2
        for dpx = 1:3
            [BIC_PERFORMANCE(b,f,dpx,c).TPR, BIC_PERFORMANCE(b,f,dpx,c).FPR] = FIND_TPRFPR(BIC_PERFORMANCE(b,f,dpx,c).bic);
        
        end
    end
end
end
LIST_TPR = zeros(2,10);
LIST_FPR = zeros(2,10);
cnt=0;
ORDER = [1 3 2 4];
for c=1:2
    for dpx=2:3
        cnt= cnt+1;
for b=1:10
    
    LIST_TPR(ORDER(cnt),b) = BIC_PERFORMANCE(b,1,dpx,c).TPR;
    LIST_FPR(ORDER(cnt),b) = BIC_PERFORMANCE(b,1,dpx,c).FPR;
end
    end
end
figure(1)
HH = plot(LIST_FPR');
set(HH(1),'Color','b','LineWidth',1,'Marker','o','LineStyle','--')
set(HH(3),'Color','b','LineWidth',1,'Marker','o','LineStyle','-')
set(HH(2),'Color','r','LineWidth',1,'Marker','+','LineStyle','--')
set(HH(4),'Color','r','LineWidth',1,'Marker','+','LineStyle','-')

legend('non-convex, T=300','convex, T=300','non-convex, T=1350','convex, T=1350','Location','northeast')
xlabel('Topology type')
ylabel('FPR')
ylim([0 0.8])
set(gca,'FontSize',11)
% saveas(figure(1),'topologyD20_diff01','epsc')
%%

function [TPR,FPR] = FIND_TPRFPR(X)
TPR = X(1)./(X(1)+X(2));
FPR = X(4)./(X(3)+X(4));
end