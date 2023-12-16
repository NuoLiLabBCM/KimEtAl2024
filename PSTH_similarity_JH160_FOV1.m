% calculating within sessions PSTH similarity
clear; clc;
nplanes = 5;
default_cd = cd;
sr = 30/nplanes;

mouse_id = 160;
FOV = 1;
case_id = 'oxoxx';
savefn = strcat('mouse',num2str(mouse_id),'_FOV',num2str(FOV),'_case_',case_id,'.mat');

session_ori{1,1} = ["2023_06_24"];  % tactile 1
session_ori{2,1} = ["2023_07_13"];  % tactile 1
session_ori{3,1} = ["2023_08_08"];  % tactile 1'
session_ori{4,1} = ["2022_12_02"];  % tactile 2'
session_ori{5,1} = ["2022_12_23"];  % auditory 1

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

for tt=1:size(sessions,1)
    disp(tt)
    cd(sessions{tt,1})
    
    cd(session_folder)
    load plane_all_matched_3sessions.mat plane_all
    cell_no = size(plane_all,1);
    
    for i=1:size(plane_all,1)
        
        % mean yes correct
        ROI_sessions{tt,1}{1,1}(i,:) = mean(plane_all{i,7},1);     % all trials
        
        n_tr = size(plane_all{i,7},1);
        % odd trials
        ROI_sessions{tt,1}{2,1}(i,:) = mean(plane_all{i,7}([1:2:n_tr],:),1);
        % even trials
        ROI_sessions{tt,1}{3,1}(i,:) = mean(plane_all{i,7}([2:2:n_tr],:),1);
        % 1st half trials
        ROI_sessions{tt,1}{4,1}(i,:) = mean(plane_all{i,7}(1:round(n_tr/2),:),1);
        % 2nd half trials
        ROI_sessions{tt,1}{5,1}(i,:) = mean(plane_all{i,7}(round(n_tr/2):end,:),1);
        
        
        % mean no correct
        ROI_sessions{tt,1}{1,2}(i,:) = mean(plane_all{i,8},1);
        
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
        ROI_sessions{tt,1}{1,3}(i,:) = mean(plane_all{i,17},1);
        
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
        ROI_sessions{tt,1}{1,4}(i,:) = mean(plane_all{i,18},1);
        
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
        ROI_sessions_ori{tt,1}{i,1}(:,:) = plane_all{i,7};
        ROI_sessions_ori{tt,1}{i,2}(:,:) = plane_all{i,8};
        ROI_sessions_ori{tt,1}{i,3}(:,:) = plane_all{i,17};
        ROI_sessions_ori{tt,1}{i,4}(:,:) = plane_all{i,18};
     
    end
    clear plane_all
    
    cd(default_cd)    
end
cd(default_cd)    


% responsiveness statistics
clear epoch    
epoch(1,:) = round([1.57 2.87 4.17 6.17]*6);
epoch = horzcat(1,epoch);    

for i=1:cell_no

    % lick left baseline activity
    for j=1:length(epoch)-1
        for tt=1:size(sessions,1)
            for jj=1:size(ROI_sessions_ori{tt,1}{i,1},1)
                stat_left{tt,j}(i,jj) = mean(ROI_sessions_ori{tt,1}{i,1}(jj,epoch(j):epoch(j+1)-1));
            end
        end
    end

    % lick right baseline activity
    for j=1:length(epoch)-1
        for tt=1:size(sessions,1)
            for jj=1:size(ROI_sessions_ori{tt,1}{i,2},1)
                stat_right{tt,j}(i,jj) = mean(ROI_sessions_ori{tt,1}{i,2}(jj,epoch(j):epoch(j+1)-1));
            end
        end
    end

    % sample-epoch selective cells
    for j=1:size(sessions,1)
        clear h p
        [h p] = ttest2(stat_left{j,2}(i,:),stat_right{j,2}(i,:));
        if p < 0.01
            stat_sample(i,j) = 1;
        else
            stat_sample(i,j) = 0;
        end
    end

    % delay-epoch selective cells
    for j=1:size(sessions,1)
        clear h p
        [h p] = ttest2(stat_left{j,3}(i,:),stat_right{j,3}(i,:));
        if p < 0.01
            stat_delay(i,j) = 1;
        else
            stat_delay(i,j) = 0;
        end
    end

    % response-epoch selective cells
    for j=1:size(sessions,1)
        clear h p
        [h p] = ttest2(stat_left{j,4}(i,:),stat_right{j,4}(i,:));
        if p < 0.01
            stat_response(i,j) = 1;
        else
            stat_response(i,j) = 0;
        end
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

stat_nonsel = horzcat(stat_nonsel_1,stat_nonsel_2,stat_nonsel_3);
stat_sel = horzcat(stat_sample,stat_delay,stat_response);

stat_nonsel_all = sum(stat_nonsel,2);
stat_sel_all = sum(stat_sel,2);


% all selective cells
    
cell_rev = find(stat_sel_all >= 0);
cell_no = length(cell_rev);

for kk=1:5
    for tt=1:size(sessions,1)
        % concatenate lick left and lick right trials
        PSTH_combined{kk,tt} = horzcat(ROI_sessions{tt,1}{kk,1}(cell_rev,:),ROI_sessions{tt,1}{kk,2}(cell_rev,:));
        % delta PSTH (lick right - lick left)
        PSTH_combined3{kk,tt} = ROI_sessions{tt,1}{kk,2}(cell_rev,:) - ROI_sessions{tt,1}{kk,1}(cell_rev,:);
    end
end   

% all non-selective cells

cell_rev2 = find(stat_sel_all == 0);
cell_no2 = length(cell_rev2);

for kk=1:5
    for tt=1:size(sessions,1)
        % concatenate lick left and lick right trials
        PSTH_combined2{kk,tt} = horzcat(ROI_sessions{tt,1}{kk,1}(cell_rev2,:),ROI_sessions{tt,1}{kk,2}(cell_rev2,:));
        % delta PSTH (lick right - lick left)
        PSTH_combined4{kk,tt} = ROI_sessions{tt,1}{kk,2}(cell_rev2,:) - ROI_sessions{tt,1}{kk,1}(cell_rev2,:);
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
        PSTH_simil_within_oe(i,jj) = r(1,2);

        clear r p
        [r p] = corrcoef(PSTH_combined{4,jj}(i,:),PSTH_combined{5,jj}(i,:));
        PSTH_simil_within_12(i,jj) = r(1,2);

    end
end

% sorting our reliable cells
for i=1:size(PSTH_simil_within_12,1)
    for j=1:size(PSTH_simil_within_12,2)
        if PSTH_simil_within_12(i,j) >= 0.5
            rel_mat(i,j) = 1;
        else
            rel_mat(i,j) = 0;
        end
    end
    
    if length(find(rel_mat(i,:) == 0)) == 0
        rel_log(i,1) = 1;
    else
        rel_log(i,1) = 0;
    end
end

rel_cell_id = find(rel_log == 1);
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




for i=1:cell_no2
    
    % similarity of concatenated PSTH (across sessionss)
    
    for kk=1
        for tt=1:size(sessions,1)
            for zz=1:size(sessions,1)
                if zz > tt
                    clear r p
                    % calculating similarity across sessionss using all trials
                    [r p] = corrcoef(PSTH_combined2{1,tt}(i,:),PSTH_combined2{1,zz}(i,:));            
                    PSTH_simil_across_all_2{tt,zz}(i,1) = r(1,2);
                    clear r p
                    % calculating similarity across sessionss using odd trials
                    [r p] = corrcoef(PSTH_combined2{2,tt}(i,:),PSTH_combined2{2,zz}(i,:));            
                    PSTH_simil_across_odd_2{tt,zz}(i,1) = r(1,2);
                    clear r p
                    % calculating similarity across sessionss using even trials
                    [r p] = corrcoef(PSTH_combined2{3,tt}(i,:),PSTH_combined2{3,zz}(i,:));            
                    PSTH_simil_across_even_2{tt,zz}(i,1) = r(1,2);
                end
            end
        end
    end
    
    for jj=1:size(sessions,1)
        clear r p
        [r p] = corrcoef(PSTH_combined2{2,jj}(i,:),PSTH_combined2{3,jj}(i,:));
        PSTH_simil_within_oe_2(i,jj) = r(1,2);

        clear r p
        [r p] = corrcoef(PSTH_combined2{4,jj}(i,:),PSTH_combined2{5,jj}(i,:));
        PSTH_simil_within_12_2(i,jj) = r(1,2);

    end
end

for tt=1:size(sessions,1)
    for zz=1:size(sessions,1)
        if zz > tt
            % selective, PSTH similarity, across sessionss, all trials
            PSTH_simil_final{1,1}(tt,zz) = mean(PSTH_simil_across_all{tt,zz}(rel_cell_id));
            % non-selective, PSTH similarity, across sessionss, all trials
            PSTH_simil_final{2,1}(tt,zz) = mean(PSTH_simil_across_all{tt,zz}(non_rel_cell_id));
            % selective, PSTH similarity, across sessionss, mean(odd vs odd, even vs even trials)
            PSTH_simil_final{1,2}(tt,zz) = (mean(PSTH_simil_across_odd{tt,zz}(rel_cell_id))+mean(PSTH_simil_across_even{tt,zz}(rel_cell_id)))/2;
            % non-selective, PSTH similarity, across sessionss, mean(odd vs odd, even vs even trials)
            PSTH_simil_final{2,2}(tt,zz) = (mean(PSTH_simil_across_odd{tt,zz}(non_rel_cell_id))+mean(PSTH_simil_across_even{tt,zz}(non_rel_cell_id)))/2;
        end
    end
end

% selecive, PSTH similarity, within sessionss (odd vs even)
PSTH_simil_final{3,1}(1,:) = mean(PSTH_simil_within_oe(rel_cell_id,:),1);
% non-selecive, PSTH similarity, within sessionss (odd vs even)
PSTH_simil_final{3,1}(2,:) = mean(PSTH_simil_within_oe(non_rel_cell_id,:),1);      

% selecive, PSTH similarity, within sessionss (1st vs 2nd half)
PSTH_simil_final{4,1}(1,:) = mean(PSTH_simil_within_12(rel_cell_id,:),1);
% non-selecive, PSTH similarity, within sessionss (1st vs 2nd half)
PSTH_simil_final{4,1}(2,:) = mean(PSTH_simil_within_12(non_rel_cell_id,:),1);        


% population selectivity vector correlation

% selective cells
% population selectivity vector correlation across sessionss
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
  
% non-selective cells
% population selectivity vector correlation across sessionss
curr = 0;
for tt=1:size(sessions,1)
    for zz=1:size(sessions,1)
        if zz > tt
            clear r p
            [r p] = corrcoef(PSTH_combined3{1,tt},PSTH_combined3{1,zz});
            pop_vec_corr_across_all_2(tt,zz) = r(1,2);
            clear r p
            [r p] = corrcoef(PSTH_combined3{2,tt}(non_rel_cell_id,:),PSTH_combined3{2,zz}(non_rel_cell_id,:));
            pop_vec_corr_across_odd_2(tt,zz) = r(1,2);
            clear r p
            [r p] = corrcoef(PSTH_combined3{3,tt}(non_rel_cell_id,:),PSTH_combined3{3,zz}(non_rel_cell_id,:));
            pop_vec_corr_across_even_2(tt,zz) = r(1,2);
            
             for yy=1:47
                 clear r p
                 [r p] = corrcoef(PSTH_combined3{1,tt}(non_rel_cell_id,yy),PSTH_combined3{1,zz}(non_rel_cell_id,yy));
                 pop_vec_corr_across_all_ts_2{tt,zz}(1,yy) = r(1,2);
                 clear r p
                 [r p] = corrcoef(PSTH_combined3{2,tt}(non_rel_cell_id,yy),PSTH_combined3{2,zz}(non_rel_cell_id,yy));
                 pop_vec_corr_across_odd_ts_2{tt,zz}(1,yy) = r(1,2);
                 clear r p
                 [r p] = corrcoef(PSTH_combined3{3,tt}(non_rel_cell_id,yy),PSTH_combined3{3,zz}(non_rel_cell_id,yy));
                 pop_vec_corr_across_even_ts_2{tt,zz}(1,yy) = r(1,2);
             end
            curr = curr + 1;
            pop_vec_corr_across_ts_re2(curr,:) = (pop_vec_corr_across_odd_ts_2{tt,zz}+pop_vec_corr_across_even_ts_2{tt,zz})/2;
            pop_vec_corr_across_oe_mean_2(tt,zz) = (pop_vec_corr_across_odd_2(tt,zz)+pop_vec_corr_across_even_2(tt,zz))/2;
        end
        if zz == tt
            clear r p
            [r p] = corrcoef(PSTH_combined3{2,tt}(non_rel_cell_id,:),PSTH_combined3{3,zz}(non_rel_cell_id,:));
            pop_vec_corr_within_2(1,zz) = r(1,2);
            clear r p
            [r p] = corrcoef(PSTH_combined3{4,tt}(non_rel_cell_id,:),PSTH_combined3{5,zz}(non_rel_cell_id,:));
            pop_vec_corr_within_2(2,zz) = r(1,2);
            
            for yy=1:47
                clear r p
                [r p] = corrcoef(PSTH_combined3{2,tt}(non_rel_cell_id,yy),PSTH_combined3{3,zz}(non_rel_cell_id,yy));
                pop_vec_corr_within_od_ts_2(zz,yy) = r(1,2);
                clear r p
                [r p] = corrcoef(PSTH_combined3{4,tt}(non_rel_cell_id,yy),PSTH_combined3{5,zz}(non_rel_cell_id,yy));
                pop_vec_corr_within_12_ts_2(zz,yy) = r(1,2);
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
take_home{1,2} = horzcat(PSTH_simil_final{3,1}(2,:)',vertcat(PSTH_simil_final{2,2}(:,2:end),zeros(1,length(sessions)-1)))';

% population vector 
take_home{2,1} = horzcat(pop_vec_corr_within(1,:)',vertcat(pop_vec_corr_across_oe_mean(:,2:end),zeros(1,length(sessions)-1)))';
take_home{2,2} = horzcat(pop_vec_corr_within_2(1,:)',vertcat(pop_vec_corr_across_oe_mean_2(:,2:end),zeros(1,length(sessions)-1)))';

save(savefn,'sessions','PSTH_combined3','ROI_sessions','PSTH_simil_across_alltoge',...
    'case_id','PSTH_simil_within_12','take_home','rel_cell_id','non_rel_cell_id')
