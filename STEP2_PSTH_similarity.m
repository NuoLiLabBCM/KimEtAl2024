%%% Kim, Daie, Li 2024. Written by Jae-Hyun Kim
%%% Calculating concatenated PSTH similarity within session 
%%% Saving mat file for further data analyses
% Sorting out ROIs from Suite2P output satisying the following criteria
% Concatenated PSTH simililarity (Type 1 dF/F0, R, 1st vs. 2nd half) >= 0.5 (OR gate)
% for the purpose of sorting out 'reliable neurons'
% Significantly selective at least one epoch from one imaging session (OR gate)

clear; clc;
nplanes = 5;        % number of imaging planes
default_cd = cd;
sr = 30/nplanes;    % sampling rate (s^-1)

% Basic information of the imaging sessions
mouse_id = 160;
FOV = 8;
case_id = 'ooxxx';
savefn = strcat('mouse',num2str(mouse_id),'_FOV',num2str(FOV),'_case_',case_id,'.mat');

session_ori{1,1} = ["2023_07_07"];  % tactile 1
session_ori{2,1} = ["2023_07_25"];  % tactile 1
session_ori{3,1} = ["xx"];  % tactile 1'
session_ori{4,1} = ["xx"];  % tactile 2'
session_ori{5,1} = ["xx"];  % auditory 1

if strcmp(case_id,'ooxxx')
    % t1-t2 / Case ooxxx
    sessions{1,1} = session_ori{1,1};
    sessions{2,1} = session_ori{2,1};
end

session_folder = [];
for z=1:size(sessions,1)
    if z == 1
        session_folder = sessions{z,1};
    else
        session_folder = strcat(session_folder,'_',sessions{z,1});
    end
end

% Importing dF/F0 dataset from imaging sessions
% Separately accumulating all trials, odd/even trials, 1st/2nd half trials
for tt=1:size(sessions,1)
    disp(tt)
    cd(sessions{tt,1})
    
    cd(session_folder)
    load plane_all_matched_3sessions.mat plane_all plane_all_decon_re
    cell_no = size(plane_all,1);
    
    for i=1:size(plane_all,1)
        
        % mean yes correct
        ROI_sessions{tt,1}{1,1}(i,:) = mean(plane_all{i,7},1);              % all trials (dF/F0)
        ROI_sessions{tt,1}{6,1}(i,:) = mean(plane_all_decon_re{i,1},1);     % all trials (deconvolved by OASIS package)
        
        n_tr = size(plane_all{i,7},1);  % total trial number
        % odd trials
        ROI_sessions{tt,1}{2,1}(i,:) = mean(plane_all{i,7}([1:2:n_tr],:),1);
        % even trials
        ROI_sessions{tt,1}{3,1}(i,:) = mean(plane_all{i,7}([2:2:n_tr],:),1);
        % 1st half trials
        ROI_sessions{tt,1}{4,1}(i,:) = mean(plane_all{i,7}(1:round(n_tr/2),:),1);
        % 2nd half trials
        ROI_sessions{tt,1}{5,1}(i,:) = mean(plane_all{i,7}(round(n_tr/2):end,:),1);
        
        
        % mean no correct
        ROI_sessions{tt,1}{1,2}(i,:) = mean(plane_all{i,8},1);              % all trials (dF/F0)
        ROI_sessions{tt,1}{6,2}(i,:) = mean(plane_all_decon_re{i,2},1);     % all trials (deconvolved by OASIS package)
        
        n_tr = size(plane_all{i,8},1);
        % odd trials
        ROI_sessions{tt,1}{2,2}(i,:) = mean(plane_all{i,8}([1:2:n_tr],:),1);
        % even trials
        ROI_sessions{tt,1}{3,2}(i,:) = mean(plane_all{i,8}([2:2:n_tr],:),1);
        % 1st half trials
        ROI_sessions{tt,1}{4,2}(i,:) = mean(plane_all{i,8}(1:round(n_tr/2),:),1);
        % 2nd half trials
        ROI_sessions{tt,1}{5,2}(i,:) = mean(plane_all{i,8}(round(n_tr/2):end,:),1);
        
        % mean yes incorrect
        ROI_sessions{tt,1}{1,3}(i,:) = mean(plane_all{i,17},1);             % all trials (dF/F0)
        ROI_sessions{tt,1}{6,3}(i,:) = mean(plane_all_decon_re{i,3},1);     % all trials (deconvolved by OASIS package)
        
        n_tr = size(plane_all{i,17},1);
        % odd trials
        ROI_sessions{tt,1}{2,3}(i,:) = mean(plane_all{i,17}([1:2:n_tr],:),1);
        % even trials
        ROI_sessions{tt,1}{3,3}(i,:) = mean(plane_all{i,17}([2:2:n_tr],:),1);
        % 1st half trials
        ROI_sessions{tt,1}{4,3}(i,:) = mean(plane_all{i,17}(1:round(n_tr/2),:),1);
        % 2nd half trials
        ROI_sessions{tt,1}{5,3}(i,:) = mean(plane_all{i,17}(round(n_tr/2):end,:),1);
        
        % mean no incorrect
        ROI_sessions{tt,1}{1,4}(i,:) = mean(plane_all{i,18},1);             % all trials (dF/F0)
        ROI_sessions{tt,1}{6,4}(i,:) = mean(plane_all_decon_re{i,4},1);     % all trials (deconvolved by OASIS package)
        
        n_tr = size(plane_all{i,18},1);
        % odd trials
        ROI_sessions{tt,1}{2,4}(i,:) = mean(plane_all{i,18}([1:2:n_tr],:),1);
        % even trials
        ROI_sessions{tt,1}{3,4}(i,:) = mean(plane_all{i,18}([2:2:n_tr],:),1);
        % 1st half trials
        ROI_sessions{tt,1}{4,4}(i,:) = mean(plane_all{i,18}(1:round(n_tr/2),:),1);
        % 2nd half trials
        ROI_sessions{tt,1}{5,4}(i,:) = mean(plane_all{i,18}(round(n_tr/2):end,:),1);
        
        % original dataset
        ROI_sessions_ori{tt,1}{i,1}(:,:) = plane_all_decon_re{i,1};
        ROI_sessions_ori{tt,1}{i,2}(:,:) = plane_all_decon_re{i,2};
        ROI_sessions_ori{tt,1}{i,3}(:,:) = plane_all_decon_re{i,3};
        ROI_sessions_ori{tt,1}{i,4}(:,:) = plane_all_decon_re{i,4};
     
    end
    clear plane_all plane_all_decon_re
    
    cd(default_cd)    
end
cd(default_cd)    

% Responsiveness statistics (non-paired t-test, p<0.01)
clear epoch epoch2    
epoch(1,:) = round([1.57 3.5 4.27]*6);      % start of each epoch
epoch2(1,:) = round([2.3 4.17 5.27]*6);     % end of each epoch
epoch = horzcat(1,epoch);    
epoch2 = horzcat(9,epoch2);    
   

for i=1:cell_no

    % lick left activity (correct trials)
    for j=1:length(epoch)
        for tt=1:size(sessions,1)
            for jj=1:floor(size(ROI_sessions_ori{tt,1}{i,1},1)/2)
                % even trials (statistics)
                stat_left{tt,j}(i,jj) = mean(ROI_sessions_ori{tt,1}{i,1}(2*jj,epoch(j):epoch2(j)));
                % odd trials (selectivity index plot, cross-validated)
                stat_left2{tt,j}(i,jj) = mean(ROI_sessions_ori{tt,1}{i,1}(2*jj-1,epoch(j):epoch2(j)));
            end
        end
    end

    % lick right activity (correct trials)
    for j=1:length(epoch)
        for tt=1:size(sessions,1)
            for jj=1:floor(size(ROI_sessions_ori{tt,1}{i,2},1)/2)
                % even trials (statistics)
                stat_right{tt,j}(i,jj) = mean(ROI_sessions_ori{tt,1}{i,2}(2*jj,epoch(j):epoch2(j)));
                % odd trials (selectivity index plot, cross-validated)
                stat_right2{tt,j}(i,jj) = mean(ROI_sessions_ori{tt,1}{i,2}(2*jj-1,epoch(j):epoch2(j)));
            end
        end
    end

    % sample-epoch selective cells
    for j=1:size(sessions,1)
        clear h p
        [h p] = ttest2(stat_left{j,2}(i,:),stat_right{j,2}(i,:));
        if p < 0.01
            stat_sample(i,j) = 1;   % significant neuron
        else
            stat_sample(i,j) = 0;
        end
        stat_sample2(i,j) = p;
    end

    % delay-epoch selective cells
    for j=1:size(sessions,1)
        clear h p
        [h p] = ttest2(stat_left{j,3}(i,:),stat_right{j,3}(i,:));
        if p < 0.01
            stat_delay(i,j) = 1;    % significant neuron
        else
            stat_delay(i,j) = 0;
        end
        stat_delay2(i,j) = p;
    end

    % response-epoch selective cells
    for j=1:size(sessions,1)
        clear h p
        [h p] = ttest2(stat_left{j,4}(i,:),stat_right{j,4}(i,:));
        if p < 0.01
            stat_response(i,j) = 1; % significant neuron
        else
            stat_response(i,j) = 0;
        end
        stat_response2(i,j) = p;
    end

    % responsive but non-selective cells
    for j=1:size(sessions,1)
        clear h p
        [h p] = ttest2(horzcat(stat_left{j,1}(i,:),stat_right{j,1}(i,:)),...
            horzcat(stat_left{j,2}(i,:),stat_right{j,2}(i,:)));
        if p < 0.01
            stat_nonsel_1(i,j) = 1;     % sample responsive, non-selective
        else
            stat_nonsel_1(i,j) = 0;
        end

        clear h p
        [h p] = ttest2(horzcat(stat_left{j,1}(i,:),stat_right{j,1}(i,:)),...
            horzcat(stat_left{j,3}(i,:),stat_right{j,3}(i,:)));
        if p < 0.01
            stat_nonsel_2(i,j) = 1;     % delay responsive, non-selective
        else
            stat_nonsel_2(i,j) = 0;
        end

        clear h p
        [h p] = ttest2(horzcat(stat_left{j,1}(i,:),stat_right{j,1}(i,:)),...
            horzcat(stat_left{j,4}(i,:),stat_right{j,4}(i,:)));
        if p < 0.01
            stat_nonsel_3(i,j) = 1;     % response responsive, non-selective
        else
            stat_nonsel_3(i,j) = 0;
        end
    end

end

% data accumulation for selectivit index plot in an epoch-specific manner
for tt=1:size(stat_left,1)
    for ttt=1:size(stat_left,2)
        resp_all{tt,ttt}(:,1) = mean(stat_left2{tt,ttt},2);
        resp_all{tt,ttt}(:,2) = mean(stat_right2{tt,ttt},2);
    end
end

stat_nonsel = horzcat(stat_nonsel_1,stat_nonsel_2,stat_nonsel_3);
stat_sel = horzcat(stat_sample,stat_delay,stat_response);
stat_p = horzcat(stat_sample2,stat_delay2,stat_response2);

stat_nonsel_all = sum(stat_nonsel,2);   % non-selective but responsive neurons (if value > 0)
stat_sel_all = sum(stat_sel,2);         % selective neurons (if value > 0)

%%% Criteria 2: Significantly selective at least one epoch from one imaging session
cell_rev = find(stat_sel_all >= 0);     % neuron ID
cell_no = length(cell_rev);             % number of neurons

for kk=1:6
    for tt=1:size(sessions,1)
        % concatenate lick left and lick right trials
        PSTH_combined{kk,tt} = horzcat(ROI_sessions{tt,1}{kk,1}(cell_rev,:),ROI_sessions{tt,1}{kk,2}(cell_rev,:));
        % delta PSTH (lick right - lick left)
        PSTH_combined3{kk,tt} = ROI_sessions{tt,1}{kk,2}(cell_rev,:) - ROI_sessions{tt,1}{kk,1}(cell_rev,:);
    end
end   


for i=1:cell_no
    
    % similarity of concatenated PSTH (across sessionss)
    
    for kk=1
        for tt=1:size(sessions,1)
            for zz=1:size(sessions,1)
                if zz > tt
                    clear r p
                    % calculating similarity across sessionss, all trials
                    [r p] = corrcoef(PSTH_combined{1,tt}(i,:),PSTH_combined{1,zz}(i,:));            
                    PSTH_simil_across_all{tt,zz}(i,1) = r(1,2);
                    clear r p
                    % calculating similarity across sessionss, odd trials
                    [r p] = corrcoef(PSTH_combined{2,tt}(i,:),PSTH_combined{2,zz}(i,:));            
                    PSTH_simil_across_odd{tt,zz}(i,1) = r(1,2);
                    clear r p
                    % calculating similarity across sessionss, even trials
                    [r p] = corrcoef(PSTH_combined{3,tt}(i,:),PSTH_combined{3,zz}(i,:));            
                    PSTH_simil_across_even{tt,zz}(i,1) = r(1,2);
                end
            end
        end
    end
    
    for jj=1:size(sessions,1)
        clear r p
        [r p] = corrcoef(PSTH_combined{2,jj}(i,:),PSTH_combined{3,jj}(i,:));
        PSTH_simil_within_oe(i,jj) = r(1,2);    % odd vs. even

        clear r p
        [r p] = corrcoef(PSTH_combined{4,jj}(i,:),PSTH_combined{5,jj}(i,:));
        PSTH_simil_within_12(i,jj) = r(1,2);    % 1st vs. 2nd half

    end
end

% sorting our reliable cells
for i=1:size(PSTH_simil_within_12,1)
    for j=1:size(PSTH_simil_within_12,2)
        if PSTH_simil_within_12(i,j) >= 0.5
            rel_mat(i,j) = 1;   % R >= 0.5
        else
            rel_mat(i,j) = 0;   % R < 0.5
        end
    end
    
    if length(find(rel_mat(i,:) == 1)) > 0
        rel_log(i,1) = 1;       % R >= 0.5 (OR gate)
    else
        rel_log(i,1) = 0;
    end
end


rel_cell_id = find(rel_log == 1);   % number of reliable neurons
non_rel_cell_id = find(rel_log == 0);

curcur = 0;
PSTH_simil_across_alltoge=[];
for tt=1:size(sessions,1)
    for zz=1:size(sessions,1)
        if zz > tt
            curcur = curcur+1;
            clear temptemp
            temptemp = horzcat(PSTH_simil_across_even{tt,zz},PSTH_simil_across_odd{tt,zz});
            PSTH_simil_across_alltoge(:,curcur) = mean(temptemp,2);
        end
    end
end

for tt=1:size(sessions,1)
    for zz=1:size(sessions,1)
        if zz > tt
            % selective, PSTH similarity, across sessionss, all trials
            PSTH_simil_final{1,1}(tt,zz) = mean(PSTH_simil_across_all{tt,zz}(rel_cell_id));
            % selective, PSTH similarity, across sessionss, mean(odd vs odd, even vs even trials)
            PSTH_simil_final{1,2}(tt,zz) = (mean(PSTH_simil_across_odd{tt,zz}(rel_cell_id))+mean(PSTH_simil_across_even{tt,zz}(rel_cell_id)))/2;
        end
    end
end

% selecive, PSTH similarity, within sessionss (odd vs even)
PSTH_simil_final{3,1}(1,:) = mean(PSTH_simil_within_oe(rel_cell_id,:),1);

% selecive, PSTH similarity, within sessionss (1st vs 2nd half)
PSTH_simil_final{4,1}(1,:) = mean(PSTH_simil_within_12(rel_cell_id,:),1);

% Population selectivity vector correlation across sessions
% similar analysis done by Schoonover et al 2021, Marks & Goard 2021
curr = 0;
for tt=1:size(sessions,1)
    for zz=1:size(sessions,1)
        if zz > tt
            clear r p
            [r p] = corrcoef(PSTH_combined3{1,tt}(rel_cell_id,:),PSTH_combined3{1,zz}(rel_cell_id,:));
            pop_vec_corr_across_all(tt,zz) = r(1,2);
            clear r p
            [r p] = corrcoef(PSTH_combined3{2,tt}(rel_cell_id,:),PSTH_combined3{2,zz}(rel_cell_id,:));
            pop_vec_corr_across_odd(tt,zz) = r(1,2);
            clear r p
            [r p] = corrcoef(PSTH_combined3{3,tt}(rel_cell_id,:),PSTH_combined3{3,zz}(rel_cell_id,:));
            pop_vec_corr_across_even(tt,zz) = r(1,2);
            
             for yy=1:47
                 clear r p
                 [r p] = corrcoef(PSTH_combined3{1,tt}(rel_cell_id,yy),PSTH_combined3{1,zz}(rel_cell_id,yy));
                 pop_vec_corr_across_all_ts{tt,zz}(1,yy) = r(1,2);
                 clear r p
                 [r p] = corrcoef(PSTH_combined3{2,tt}(rel_cell_id,yy),PSTH_combined3{2,zz}(rel_cell_id,yy));
                 pop_vec_corr_across_odd_ts{tt,zz}(1,yy) = r(1,2);
                 clear r p
                 [r p] = corrcoef(PSTH_combined3{3,tt}(rel_cell_id,yy),PSTH_combined3{3,zz}(rel_cell_id,yy));
                 pop_vec_corr_across_even_ts{tt,zz}(1,yy) = r(1,2);
             end
            curr = curr + 1;
            pop_vec_corr_across_ts_re(curr,:) = (pop_vec_corr_across_odd_ts{tt,zz}+pop_vec_corr_across_even_ts{tt,zz})/2;             
            pop_vec_corr_across_oe_mean(tt,zz) = (pop_vec_corr_across_odd(tt,zz)+pop_vec_corr_across_even(tt,zz))/2;
        end
        if zz == tt
            clear r p
            [r p] = corrcoef(PSTH_combined3{2,tt}(rel_cell_id,:),PSTH_combined3{3,zz}(rel_cell_id,:));
            pop_vec_corr_within(1,zz) = r(1,2);
            clear r p
            [r p] = corrcoef(PSTH_combined3{4,tt}(rel_cell_id,:),PSTH_combined3{5,zz}(rel_cell_id,:));
            pop_vec_corr_within(2,zz) = r(1,2);
            
            for yy=1:47
                clear r p
                [r p] = corrcoef(PSTH_combined3{2,tt}(rel_cell_id,yy),PSTH_combined3{3,zz}(rel_cell_id,yy));
                pop_vec_corr_within_od_ts(zz,yy) = r(1,2);
                clear r p
                [r p] = corrcoef(PSTH_combined3{4,tt}(rel_cell_id,yy),PSTH_combined3{5,zz}(rel_cell_id,yy));
                pop_vec_corr_within_12_ts(zz,yy) = r(1,2);
            end
        end
    end
end  

disp(strcat('total=',num2str(length(rel_log))))
disp(strcat('reliable=',num2str(length(rel_cell_id))))
disp(strcat('non-reliable=',num2str(length(non_rel_cell_id))))

% final data reorganization

% PSTH similarity vector
take_home{1,1} = horzcat(PSTH_simil_final{3,1}(1,:)',vertcat(PSTH_simil_final{1,2}(:,2:end),zeros(1,length(sessions)-1)))';

% population vector 
take_home{2,1} = horzcat(pop_vec_corr_within(1,:)',vertcat(pop_vec_corr_across_oe_mean(:,2:end),zeros(1,length(sessions)-1)))';

save(savefn,'sessions','PSTH_combined3','ROI_sessions','PSTH_simil_across_alltoge',...
    'case_id','PSTH_simil_within_12','take_home','rel_cell_id','non_rel_cell_id','stat_sel','stat_nonsel','stat_p','resp_all')
