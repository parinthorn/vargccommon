% This experiment calculate ROC curve of estimated model both convex and non-convex.

clear

SFolderName = 'G:\Shared drives\MASTER_DRIVE\THESIS\RESULTS\FormulationC_Model16Oct19_d20';
DataFolderName = 'G:\Shared drives\MASTER_DRIVE\THESIS\THESIS_DATA_conference';
ModelFolderName = 'G:\Shared drives\MASTER_DRIVE\THESIS\MODEL';
Resultd10='G:\Shared drives\MASTER_DRIVE\THESIS\RESULTS\FormulationC_Model16Oct19_d10';
Resultd20='G:\Shared drives\MASTER_DRIVE\THESIS\RESULTS\FormulationC_Model16Oct19_d20';
Modeld10= load(strcat(ModelFolderName,'\MODEL16Oct19_d10'));

Modeld20= load(strcat(ModelFolderName,'\MODEL16Oct19_d20'));
% Result.d10 = struct();
% Result.d20 = struct();
DP = [50,300,1350];
diff_density = [0.01,0.05];
for f=1:2
    for dpx=1:3
        ROC1 = [0 0 0 0];
        ROC2 = [0 0 0 0];
                ROC3 = [0 0 0 0];
        ROC4 = [0 0 0 0];
        for b=1:10
            load(strcat(Resultd10,'\ESTIMATED_BANK_',int2str(b),'_D',int2str(f),'_',int2str(DP(dpx))))
            Ed10=E;
            
            load(strcat(Resultd20,'\ESTIMATED_BANK_',int2str(b),'_D',int2str(f),'_',int2str(DP(dpx))))
            Ed20=E;
            load(strcat(Resultd10,'\ESTIMATED_BANK_',int2str(b),'_D',int2str(f),'_',int2str(DP(dpx)),'_CVX'))
            Ed10cvx=E;
            
            load(strcat(Resultd20,'\ESTIMATED_BANK_',int2str(b),'_D',int2str(f),'_',int2str(DP(dpx)),'_CVX'))
            Ed20cvx=E;
            for ii=1:length(E)
                O1 = Ed10(ii).M.stat.consistency;
                O2 = Ed20(ii).M.stat.consistency;
                O3 = Ed10cvx(ii).M.stat.consistency;
                O4 = Ed20cvx(ii).M.stat.consistency;                
                ROC1 = ROC1 + O1(1,:,:,:);
                ROC2 = ROC2 + O2(1,:,:,:);
                ROC3 = ROC3 + O3(1,:,:,:);
                ROC4 = ROC4 + O4(1,:,:,:);
                
            end
            
            %             Atrued10=Modeld10.E(b).M.(strcat('D',int2str(f)));
            %             Atrued20=Modeld20.E(b).M.(strcat('D',int2str(f)));
            
        end
        Result(f,dpx).ROCd10 = ROC1;
        Result(f,dpx).ROCd20 = ROC2;
        Result(f,dpx).ROCd10cvx = ROC3;
        Result(f,dpx).ROCd20cvx = ROC4;
        
        Result(f,dpx).info.Datapoint = DP(dpx);
        Result(f,dpx).info.diff_density =diff_density(f);
    end
end