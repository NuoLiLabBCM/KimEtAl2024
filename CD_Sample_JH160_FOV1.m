% cross-session CD analysis (deconvolved)
% additional criteria for data analysis
% 1. Trial number >= 50
% 2. Inferred spike rate <= 15 (every bin)
% 3. CD window: 1 sec (6 bins)

clear;clc;
nplanes = 5;
addpath('D:\Jae-Hyun Kim\Ca_imaging_functions_JH')
addpath('D:\Jae-Hyun Kim\CD analysis\func')

session_ori{1,1} = ["2023_06_24"];  % tactile 1
session_ori{2,1} = ["2023_07_13"];  % tactile 1
session_ori{3,1} = ["2023_08_08"];  % tactile 1'
session_ori{4,1} = ["2022_12_02"];  % tactile 2'
session_ori{5,1} = ["2022_12_23"];  % auditory 1

mouse_id = 160;
FOV = 1;
case_id = 'oxoxx';
rep_no = 3;
gatetype = 1;   %1 for OR, 2 for AND

if gatetype == 1
    savefnn = strcat('mouse',num2str(mouse_id),'_FOV',num2str(FOV),'_case_',case_id,'_OR_sample_',num2str(rep_no),'.mat');
    savefign1 = strcat('figure_mouse',num2str(mouse_id),'_FOV',num2str(FOV),'_case_',case_id,'_OR_sample_single_trial.fig');
    savefign2 = strcat('figure_mouse',num2str(mouse_id),'_FOV',num2str(FOV),'_case_',case_id,'_OR_sample_CDaxis_',num2str(rep_no),'.fig');
elseif gatetype == 2
    savefnn = strcat('mouse',num2str(mouse_id),'_FOV',num2str(FOV),'_case_',case_id,'_AND_sample_',num2str(rep_no),'.mat');
    savefign1 = strcat('figure_mouse',num2str(mouse_id),'_FOV',num2str(FOV),'_case_',case_id,'_AND_sample_single_trial.fig');
    savefign2 = strcat('figure_mouse',num2str(mouse_id),'_FOV',num2str(FOV),'_case_',case_id,'_AND_sample_CDaxis_',num2str(rep_no),'.fig');
end

casefn = strcat('mouse',num2str(mouse_id),'_FOV',num2str(FOV),'_case_',case_id,'.mat');
load(casefn, 'PSTH_simil_within_12')
% identify reliable cells (AND gate or ORgate)

for i=1:size(PSTH_simil_within_12,1)
    for j=1:size(PSTH_simil_within_12,2)
        if PSTH_simil_within_12(i,j) >= 0.5
            rel_mat(i,j) = 1;
        else
            rel_mat(i,j) = 0;
        end
    end
    
    if gatetype == 1    % OR gate
        
        if isempty(find(rel_mat(i,:) == 1))
            rel_log(i,1) = 0;
        else
            rel_log(i,1) = 1;
        end
        
    elseif gatetype == 2    % AND gate
        
        if isempty(find(rel_mat(i,:) == 0))
            rel_log(i,1) = 1;
        else
            rel_log(i,1) = 0;
        end
    end
end

real_id = find(rel_log == 1);
cell_no = length(real_id);

if strcmp(case_id,'ooxxx')
    % t1-t2 / Case ooxxx
    sessions{1,1} = session_ori{1,1};
    sessions{2,1} = session_ori{2,1};
    
elseif strcmp(case_id,'oooxx')
    % t1-t2-t1' / Case oooxx
    sessions{1,1} = session_ori{1,1};
    sessions{2,1} = session_ori{2,1};
    sessions{3,1} = session_ori{3,1};
    
elseif strcmp(case_id,'oxoxx')
    % t1-t1' / Case oxoxx
    sessions{1,1} = session_ori{1,1};
    sessions{2,1} = session_ori{3,1};
    
elseif strcmp(case_id,'oxxox')
    % t1-t2' / Case oxxox
    sessions{1,1} = session_ori{1,1};
    sessions{2,1} = session_ori{4,1};
    
elseif strcmp(case_id,'xooxx')
    % t2-t1' / Case xooxx
    sessions{1,1} = session_ori{2,1};
    sessions{2,1} = session_ori{3,1};
    
elseif strcmp(case_id,'xooox')
    % t2-t1'-t2' / Case xooox
    sessions{1,1} = session_ori{2,1};
    sessions{2,1} = session_ori{3,1};
    sessions{3,1} = session_ori{4,1};
    
elseif strcmp(case_id,'xxoox')
    % t1'-t2' / Case xxoox
    sessions{1,1} = session_ori{3,1};
    sessions{2,1} = session_ori{4,1};
    
elseif strcmp(case_id,'ooxxx')
    % t2-t2' / Case xoxox
    sessions{1,1} = session_ori{2,1};
    sessions{2,1} = session_ori{4,1};
    
elseif strcmp(case_id,'oooox')
    % t1-t2-t1'-t2' / Case oooox
    sessions{1,1} = session_ori{1,1};
    sessions{2,1} = session_ori{2,1};
    sessions{3,1} = session_ori{3,1};
    sessions{4,1} = session_ori{4,1};
    
elseif strcmp(case_id,'xxooo')
    % t2-t1'-t2'-a / Case xxooo
    sessions{1,1} = session_ori{3,1};
    sessions{2,1} = session_ori{4,1};
    sessions{3,1} = session_ori{5,1};
    
elseif strcmp(case_id,'oxxxo')
    % t1-a / Case oxxxo
    sessions{1,1} = session_ori{1,1};
    sessions{2,1} = session_ori{5,1};
    
elseif strcmp(case_id,'xoxxo')
    % t2-a / Case xoxxo
    sessions{1,1} = session_ori{2,1};
    sessions{2,1} = session_ori{5,1};
    
elseif strcmp(case_id,'xxoxo')
    % t1'-a / Case xxoxo
    sessions{1,1} = session_ori{3,1};
    sessions{2,1} = session_ori{5,1};
    
elseif strcmp(case_id,'xxxoo')
    % t2'-a / Case xxxoo
    sessions{1,1} = session_ori{4,1};
    sessions{2,1} = session_ori{5,1};
end

session_folder = [];
for z=1:size(sessions,1)
    if z == 1
        session_folder = sessions{z,1};
    else
        session_folder = strcat(session_folder,'_',sessions{z,1});
    end
end

% load datasets from individual session

default_cd = cd;
sr = 30/nplanes;
for tt=1:size(sessions,1)
    
    % move to specific session
    cd(sessions{tt,1}) 
    
    cd(session_folder)
    
    load plane_all_matched_3sessions.mat plane_all_decon_re 
    
    % dataset reorganization (neuron x time x trial)
    for zz=1:length(real_id)
        dataset_ori{tt,1}(zz,:,:) = plane_all_decon_re{real_id(zz),1}';
        dataset_ori{tt,2}(zz,:,:) = plane_all_decon_re{real_id(zz),2}';
        dataset_ori{tt,3}(zz,:,:) = plane_all_decon_re{real_id(zz),3}';
        dataset_ori{tt,4}(zz,:,:) = plane_all_decon_re{real_id(zz),4}';
    end
    
    clear plane_all
    cd(default_cd)
end

% convert NaN to zero (<=15)
for tt=1:size(dataset_ori,1)
    dataset_ori{tt,1}(isnan(dataset_ori{tt,1})) = 0;
    dataset_ori{tt,2}(isnan(dataset_ori{tt,2})) = 0;
        
    %dataset_ori{tt,1}(find(dataset_ori{tt,1} > 15)) = 0;
    %dataset_ori{tt,2}(find(dataset_ori{tt,2} > 15)) = 0;
end

% trial number (>=50)
for tt=1:size(dataset_ori,1)
    para_trial_no(tt,1) = size(dataset_ori{tt,1},3);  
    disp(strcat(num2str(para_trial_no(tt,1)),'(left)-',sessions{tt,1}));
    
    para_trial_no(tt,2) = size(dataset_ori{tt,2},3);
    disp(strcat(num2str(para_trial_no(tt,2)),'(right)-',sessions{tt,1}));
    
    % 0 possible, 1 discard
    para_trial_no_cri(1,tt) = length(find(para_trial_no(tt,:) < 50));
end

% inferred spike rate (0< & <=15)
for tt=1:size(dataset_ori,1)
    for ttt=1:length(real_id)
        clear temp_t1 temp_t2 temp_t3 temp_t4
        temp_t1(:,:) = dataset_ori{tt,1}(ttt,:,:);
        temp_t2(:,:) = dataset_ori{tt,2}(ttt,:,:);
        temp_t3 = mean(temp_t1,1);
        temp_t4 = mean(temp_t2,1);
        
        para_inferred_spk{tt,1}(ttt,1) = max(temp_t1,[],'all');
        para_inferred_spk{tt,1}(ttt,2) = max(temp_t2,[],'all');
%         para_inferred_spk{tt,2}(ttt,1) = 0;
%         
%         if max(temp_t3) == temp_t3(1)
%             para_inferred_spk{tt,2}(ttt,1) = 1;
%         end
%         if max(temp_t4) == temp_t4(1)
%             para_inferred_spk{tt,2}(ttt,1) = 1;
%         end
        
        % 0 include, >0 discard
        %para_inferred_spk_cri(ttt,tt) = length(find(para_inferred_spk{tt,1}(ttt,:)>15)) + length(find(para_inferred_spk{tt,1}(ttt,:)==0)) + para_inferred_spk{tt,2}(ttt,1);
        para_inferred_spk_cri(ttt,tt) = length(find(para_inferred_spk{tt,1}(ttt,:)>15)) + length(find(para_inferred_spk{tt,1}(ttt,:)==0));
                
    end
end


para_inferred_spk_cri = sum(para_inferred_spk_cri,2);
if length(find(para_inferred_spk_cri > 0)) == 0
    disp('no weird inferred spk rate')
elseif length(find(para_inferred_spk_cri > 0)) > 0    
    disp('weird inferred spk rate detected')
    weird_units = find(para_inferred_spk_cri > 0);
    
    dataset_ori_temp = dataset_ori;
    clear dataset_ori

    curr = 0;
    for tttt=1:length(real_id)
        if sum(ismember(weird_units,tttt)) == 0  % normal units
            curr = curr + 1;
            for tt=1:size(dataset_ori_temp,1)
                for ttt=1:size(dataset_ori_temp,2)
                    dataset_ori{tt,ttt}(curr,:,:) = dataset_ori_temp{tt,ttt}(tttt,:,:);
                end
            end
        elseif sum(ismember(weird_units,tttt)) > 0  % weird units
            disp(strcat('weird unit:',num2str(tttt)))                
        end
    end
    
    clear dataset_ori_temp
    cell_no = curr;
end
    


total_timebin = size(dataset_ori{tt,1},2);
time_bins = 1:1:total_timebin;

% define CD time window (CD sample, CD delay, CD response)
% 5 bins
epoch_bin = floor([1.57 3.17 4.37]*sr);    % use the same time period (1.5s)
i_timebin_s = [epoch_bin(1):1:epoch_bin(1)+5];    % sample CD time window
i_timebin_d = [epoch_bin(2):1:epoch_bin(2)+5];    % delay CD time window
i_timebin_r = [epoch_bin(3):1:epoch_bin(3)+8];    % response CD time window

% iteration of train and test
% train and test trials should be non-overlapped

for tt=1:size(sessions,1)
    
    disp(num2str(tt))
    disp(strcat('cell no_',num2str(cell_no)))
    tic
    
    % trial info sorting out
    for ttt=1:size(sessions,1)  % test session
        
        % tt == ttt: decoder within session
        % tt ~= ttt: decoder across sessions
        
        % interation of train and test
        for zzt = 1:rep_no
            
            % trial no. (train)
            yes_n = size(dataset_ori{tt,2},3);
            no_n = size(dataset_ori{tt,1},3);
            
            % trial no. (test)
            yes_n2 = size(dataset_ori{ttt,2},3);
            no_n2  = size(dataset_ori{ttt,1},3);
            
            % train dataset trial info
            train_trial_info{tt,1}{zzt,1} = randsample([1:1:yes_n],round(yes_n/2));
            train_trial_info{tt,1}{zzt,2} = randsample([1:1:no_n],round(no_n/2));
            
            % test dataset trial info
            if tt == ttt   % within session: sorting out non-overlapping trials                
                test_trial_info{tt,ttt}{zzt,1} = find(ismember([1:1:yes_n],train_trial_info{tt,1}{zzt,1}) == 0);
                test_trial_info{tt,ttt}{zzt,2} = find(ismember([1:1:no_n],train_trial_info{tt,1}{zzt,2}) == 0);
                
            elseif tt ~= ttt   % across session: simply sorting out half trials
                test_trial_info{tt,ttt}{zzt,1} = randsample([1:1:yes_n2],round(yes_n2/2));
                test_trial_info{tt,ttt}{zzt,2} = randsample([1:1:no_n2],round(no_n2/2));
            end
        end
    end
    
    
    % train decoder
    
    for zzt = 1:rep_no
            
        i_yes_correct_train = train_trial_info{tt,1}{zzt,1};
        i_no_correct_train = train_trial_info{tt,1}{zzt,2};

        for timebin = 1:total_timebin

            % training data [trial x neuron]
            data_yes = squeeze(dataset_ori{tt,2}(:,timebin,i_yes_correct_train))';
            data_no = squeeze(dataset_ori{tt,1}(:,timebin,i_no_correct_train))';

            % LDA calculation
            w_iTimeBin = KD_LDA2(data_yes,data_no);
            % first LDA axis
            w_iTimeBin = w_iTimeBin(:,1);
            % accumulating first LDA axis per time bin
            w_LDA_correct{tt,zzt}(:,timebin) = w_iTimeBin;
        end

        % average LDA (epoch-specific)
        w_LDA_mean_s{1,zzt} = mean(w_LDA_correct{tt,zzt}(:,i_timebin_s)')';

        % variance in activity (epoch-specific, trial-averaged)
        Rt_co_s{tt,zzt} = mean(dataset_ori{tt,2}(:,i_timebin_s,i_yes_correct_train),3);
        Lt_co_s{tt,zzt} = mean(dataset_ori{tt,1}(:,i_timebin_s,i_no_correct_train),3);

        activityRL_s{tt,zzt} = [Lt_co_s{tt,zzt} Rt_co_s{tt,zzt}];

        clear u_s s_s v_s
        
        [u_s s_s v_s] = svd(activityRL_s{tt,1}');

        % selRL' = u * s * v';
        % u: left singular vector [time bin x time bin]
        % s: singular vector [time bin x Neuron] / diagonal matrix
        % v: right singular vector  [Neuron x Neuron]

        % transformation vector
        orthonormal_basis_s{tt,zzt} = Gram_Schmidt_process([w_LDA_mean_s{1,zzt} v_s]);

        proj_allDim_s{tt,zzt} = activityRL_s{tt,zzt}'*orthonormal_basis_s{tt,zzt};

        var_allDim_s{tt,zzt} = sum(proj_allDim_s{tt,zzt}.^2);

    end
    
    
    % test decoder accuracy
    
    for ttt=1:size(sessions,1)
        
        for zzt = 1:rep_no
            
            yes_correct_test_trial = test_trial_info{tt,ttt}{zzt,1};            
            no_correct_test_trial = test_trial_info{tt,ttt}{zzt,2};
            
            yes_correct_train_trial = train_trial_info{tt,1}{zzt,1};            
            no_correct_train_trial = train_trial_info{tt,1}{zzt,2};
            
            % trial averaged (test)
            yes_temp = mean(dataset_ori{ttt,2}(:,:,yes_correct_test_trial),3);
            no_temp = mean(dataset_ori{ttt,1}(:,:,no_correct_test_trial),3);
            
            % single trial dataset (test)
            yes_temp_trial = dataset_ori{ttt,2}(:,:,yes_correct_test_trial);
            no_temp_trial = dataset_ori{ttt,1}(:,:,no_correct_test_trial);
        
            % projection of trial-averaged trace (test)
            CD_proj_s_test{tt,ttt}{zzt,1} = yes_temp'*orthonormal_basis_s{tt,zzt};
            CD_proj_s_test{tt,ttt}{zzt,2} = no_temp'*orthonormal_basis_s{tt,zzt};
            
            % projection of single trial dataset (test)
            
            for zzz=1:size(yes_temp_trial,3)
                CD_proj_s_trial_test{tt,ttt}{zzt,1}{1,1}(:,:,zzz) = squeeze(yes_temp_trial(:,:,zzz))'*orthonormal_basis_s{tt,zzt};
                CD_proj_s_trial_test{tt,ttt}{zzt,1}{2,1}(zzz,:) = squeeze(CD_proj_s_trial_test{tt,ttt}{zzt,1}{1,1}(:,1,zzz));
            end

            for zzz=1:size(no_temp_trial,3)
                CD_proj_s_trial_test{tt,ttt}{zzt,1}{1,2}(:,:,zzz) = squeeze(no_temp_trial(:,:,zzz))'*orthonormal_basis_s{tt,zzt};
                CD_proj_s_trial_test{tt,ttt}{zzt,1}{2,2}(zzz,:) = squeeze(CD_proj_s_trial_test{tt,ttt}{zzt,1}{1,2}(:,1,zzz));
            end       
            
            % trial averaged (train)
            yes_temp = mean(dataset_ori{tt,2}(:,:,yes_correct_train_trial),3);
            no_temp = mean(dataset_ori{tt,1}(:,:,no_correct_train_trial),3);
            
            % single trial dataset (train)
            yes_temp_trial = dataset_ori{tt,2}(:,:,yes_correct_train_trial);
            no_temp_trial = dataset_ori{tt,1}(:,:,no_correct_train_trial);
        
            % projection of trial-averaged trace (train)
            CD_proj_s_train{tt,ttt}{zzt,1} = yes_temp'*orthonormal_basis_s{tt,zzt};
            CD_proj_s_train{tt,ttt}{zzt,2} = no_temp'*orthonormal_basis_s{tt,zzt};
            
            % projection of single trial dataset (train)
            
            for zzz=1:size(yes_temp_trial,3)
                CD_proj_s_trial_train{tt,ttt}{zzt,1}{1,1}(:,:,zzz) = squeeze(yes_temp_trial(:,:,zzz))'*orthonormal_basis_s{tt,zzt};
                CD_proj_s_trial_train{tt,ttt}{zzt,1}{2,1}(zzz,:) = squeeze(CD_proj_s_trial_train{tt,ttt}{zzt,1}{1,1}(:,1,zzz));
            end

            for zzz=1:size(no_temp_trial,3)
                CD_proj_s_trial_train{tt,ttt}{zzt,1}{1,2}(:,:,zzz) = squeeze(no_temp_trial(:,:,zzz))'*orthonormal_basis_s{tt,zzt};
                CD_proj_s_trial_train{tt,ttt}{zzt,1}{2,2}(zzz,:) = squeeze(CD_proj_s_trial_train{tt,ttt}{zzt,1}{1,2}(:,1,zzz));
            end 
            
            
        end
        
    end
    toc   
    
end

i_mode = 1;

time_bins_ticks = [1:1:length(time_bins)]/sr - 4.17;
epoch = [1.57 2.87 4.17]-4.17;

cd('Figures')
figure
set(gcf,'position',[650 50 600 700])
no_se = size(sessions,1);
for tt=1:no_se
    for zz=1:no_se
        subplot(no_se,no_se,no_se*(tt-1)+zz)
        hold on
        plot(time_bins_ticks, CD_proj_s_trial_test{tt,zz}{1,1}{2,1},'color',[0 0 1 .2],'LineWidth',.1);
        hold on
        plot(time_bins_ticks, CD_proj_s_trial_test{tt,zz}{1,1}{2,2},'color',[1 0 0 .2],'LineWidth',.1);
        hold on
        plot(time_bins_ticks, CD_proj_s_test{tt,zz}{1,1}(:,i_mode),'color','b','LineWidth',1);
        hold on
        plot(time_bins_ticks, CD_proj_s_test{tt,zz}{1,2}(:,i_mode),'color','r','LineWidth',1);
        
        maxlim = max(vertcat(CD_proj_s_test{tt,zz}{1,1}(:,i_mode), CD_proj_s_test{tt,zz}{1,2}(:,i_mode)));
        minlim = min(vertcat(CD_proj_s_test{tt,zz}{1,1}(:,i_mode), CD_proj_s_test{tt,zz}{1,2}(:,i_mode)));
        xlim([min(time_bins_ticks) max(time_bins_ticks)])
        %ylim([minlim maxlim])
        for i=1:3
            line([epoch(i) epoch(i)], [minlim maxlim],'color','k')
        end
        xlabel('time (s)')
        ylabel('CD projection (a.u.)')
        hold on
        title(strcat(sessions{tt,1},'=>',sessions{zz,1}),'interpreter','none')
    end
end

hold on
sgtitle('Sample CD - deconvolved')

savefig(savefign1)

shuffle_rep = 1;
no_se = size(sessions,1);
for zzt = 1:rep_no
    for tt=1:no_se
        for zz=1:no_se
            % train
            CD_proj_s_tr_sct_train{tt,zz}{zzt,1} = mean(CD_proj_s_trial_train{tt,zz}{zzt,1}{2,1}(:,i_timebin_s),2);
            CD_proj_s_tr_sct_train{tt,zz}{zzt,2} = mean(CD_proj_s_trial_train{tt,zz}{zzt,1}{2,2}(:,i_timebin_s),2);
            
            % test
            CD_proj_s_tr_sct_test{tt,zz}{zzt,1} = mean(CD_proj_s_trial_test{tt,zz}{zzt,1}{2,1}(:,i_timebin_s),2);
            CD_proj_s_tr_sct_test{tt,zz}{zzt,2} = mean(CD_proj_s_trial_test{tt,zz}{zzt,1}{2,2}(:,i_timebin_s),2);
        end
        
        % mean of lick right/left correct points on CD response projection
        CD_proj_s_DB_m{zzt,1}(tt,1) = mean(CD_proj_s_tr_sct_train{tt,tt}{zzt,1});
        CD_proj_s_DB_m{zzt,1}(tt,2) = mean(CD_proj_s_tr_sct_train{tt,tt}{zzt,2});
        
        % variance of lick right/left correct points on CD response projections
        CD_proj_s_DB_v{zzt,1}(tt,1) = (std(CD_proj_s_tr_sct_train{tt,tt}{zzt,1}))^2;
        CD_proj_s_DB_v{zzt,1}(tt,2) = (std(CD_proj_s_tr_sct_train{tt,tt}{zzt,2}))^2;
        
        % final decision boundary
        CD_proj_s_DB(tt,zzt) = (CD_proj_s_DB_m{zzt,1}(tt,1)/CD_proj_s_DB_v{zzt,1}(tt,1) + CD_proj_s_DB_m{zzt,1}(tt,2)/CD_proj_s_DB_v{zzt,1}(tt,2)) / (1/CD_proj_s_DB_v{zzt,1}(tt,1) + 1/CD_proj_s_DB_v{zzt,1}(tt,2));
        
        
        % check classifier accuracy using single trials within/across sessions
        for zz=1:no_se
            
            % CD sample classifier performance using early sample dF/F0
            clear temp_a temp_b temp_c temp_d
            % lick right trials
            temp_a = CD_proj_s_tr_sct_test{tt,zz}{zzt,1} - CD_proj_s_DB(tt,zzt);
            % lick left trials
            temp_b = CD_proj_s_tr_sct_test{tt,zz}{zzt,2} - CD_proj_s_DB(tt,zzt);
            
            temp_c = vertcat(temp_a,temp_b);
            
            for yyy=1:shuffle_rep
                clear temp_d
                temp_d = randsample(temp_c,round(length(temp_c)/3));
                temp_e(yyy,1) = length(find(temp_d > 0)) / length(temp_d);
            end
            CD_proj_s_DB_classifier_s_shuffle(tt,zz,zzt) = mean(temp_e);
            clear temp_e
            
            CD_proj_s_DB_classifier_s(tt,zz,zzt) = (length(find(temp_a > 0)) + length(find(temp_b < 0))) / (length(temp_a)+length(temp_b));
            
        end
    end
end


decoder_s_s(:,:) = mean(CD_proj_s_DB_classifier_s,3);

for i=1:length(sessions)
    for t=1:rep_no
        axis_all{1,i}(t,:) = orthonormal_basis_s{i,t}(:,1);
        axis_all{2,i}(t,:) = orthonormal_basis_s{i,t}(:,2);
        axis_all{3,i}(t,:) = orthonormal_basis_s{i,t}(:,3);
    end
end


for j=1:length(sessions)
    curr2 = 0;
    for jj=1:length(sessions)   
        if j == jj  
            curr1 = 0;
            for i=1:rep_no
                for ii=1:rep_no
                    if i < ii
                        curr1 = curr1 + 1;
                        dot_all{j,jj}(curr1,1) = dot(axis_all{1,j}(i,:),axis_all{1,jj}(ii,:));
                        dot_all2{j,jj}(curr1,1) = dot(axis_all{2,j}(i,:),axis_all{2,jj}(ii,:));
                        dot_all3{j,jj}(curr1,1) = dot(axis_all{3,j}(i,:),axis_all{3,jj}(ii,:));
                    end
                end
            end
            
        elseif j ~= jj
            curr2 = 0;
            for i=1:rep_no
                for ii=1:rep_no
                    curr2 = curr2 + 1;
                    dot_all{j,jj}(curr2,1) = dot(axis_all{1,j}(i,:),axis_all{1,jj}(ii,:));
                    dot_all2{j,jj}(curr2,1) = dot(axis_all{2,j}(i,:),axis_all{2,jj}(ii,:));
                    dot_all3{j,jj}(curr2,1) = dot(axis_all{3,j}(i,:),axis_all{3,jj}(ii,:));
                end
                    
            end
        end
    end
end

figure()
tt=0;
sz = 3;
for j=1:length(sessions)
    for jj=1:length(sessions)
        tt = tt+1;
        subplot(length(sessions),length(sessions),tt)
        scatter(rand(1,length(dot_all{j,jj})),dot_all{j,jj},sz,'k')
        hold on
        %scatter(rand(1,length(dot_all2{j,jj})),dot_all2{j,jj},sz,'r')
        %hold on
        %scatter(rand(1,length(dot_all3{j,jj})),dot_all3{j,jj},sz,'b')
        %hold on
        title(mean(dot_all{j,jj}))
        ylim([-1 1])
    end
end
savefig(savefign2)

cd(default_cd)

save(savefnn,'cell_no','rep_no','shuffle_rep','decoder_s_s','CD_proj_s_test',...
    'dot_all','axis_all','dot_all2','dot_all3','orthonormal_basis_s')
