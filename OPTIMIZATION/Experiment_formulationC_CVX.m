%% This experiment estimate VAR with formulation C by ADMM (CONVEX CASE)
clear
clc
SFolderName = 'G:\Shared drives\MASTER_DRIVE\THESIS\RESULTS\FormulationC_Model16Oct19_d20';
DataFolderName = 'G:\Shared drives\MASTER_DRIVE\THESIS\THESIS_DATA_conference';
ModelFolderName = 'G:\Shared drives\MASTER_DRIVE\THESIS\MODEL';
Mtype = {'C','D','S'};
dp = [50,300,1350];
[P,~] = offdiagJSS(15,2,4);
for b = 1:10 %C12
    s = load(strcat(ModelFolderName,'\MODEL16Oct19_d20.mat'));
    fprintf('Databank number: %d \n',b)
    for tp = 2:2
        Mname = Mtype{tp};
        for f=1:2
            disp(strcat('Model name:',Mname,int2str(f)))
            r = load(strcat(DataFolderName,'\DATA_BANK_',int2str(b),'_',Mname,int2str(f),'_d20.mat'));
%             trials = size(r.data.y,4);
            trials = 4;
            for dpx = 1:3
                E(trials).M =struct();
                Num = dp(dpx);
                fprintf('Data point : %d \n',Num)
                parfor t=1:trials
                    M = est_Formulation_c_CVX(r.data.y(:,1:Num,:,t),s.E(b).M.(strcat(Mname,int2str(f))).A,P);
                    E(t).M = M;
                    E(t).M.DataInfo = r.data.info;
                    E(t).M.DataSize = Num;
                end
                save(strcat(SFolderName,'\ESTIMATED_BANK_',int2str(b),'_',Mname,int2str(f),'_',int2str(Num),'_CVX','.mat'),'E')
                clear E
            end
        end
    end
end