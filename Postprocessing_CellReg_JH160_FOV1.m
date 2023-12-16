clear;clc
default_cd = cd;

session_ori{1,1} = ["2023_06_24"];  % tactile 1
session_ori{2,1} = ["2023_07_13"];  % tactile 1
session_ori{3,1} = ["2023_08_08"];  % tactile 1'
session_ori{4,1} = ["2022_12_02"];  % tactile 2'
session_ori{5,1} = ["2022_12_23"];  % auditory 1

case_id = 'oxoxx';

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

cd('Multi_sessions')
cd(session_folder)

currdir = dir;
clc;
curr_cd = cd;
for i=1:5
    cd(strcat('plane',num2str(i)))
    indir = dir;
    if strcmp(indir(end).name,'single_cell_matched.mat')
        load single_cell_matched.mat
        final_cell_id{i,1} = live_roi_sessions;
        disp(strcat('plane_',num2str(i),'_',num2str(length(live_roi_sessions))))
    else
        final_cell_id{i,1} = [];
    end
    cd(curr_cd)
end

save final_cell_id.mat final_cell_id

tttt= vertcat(final_cell_id{:});
disp(strcat('combined_',num2str(size(tttt,1))))
clear ans
cd(default_cd)

session_order=[];
for tt = 1:size(sessions,1)
    if size(sessions,1) == 3
        session_order = [2 1 3];
    elseif size(sessions,1) == 2
        session_order = [1 2];
    elseif size(sessions,1) == 4
        session_order = [];
    end
    
    cd(sessions{session_order(tt),1})
    clearvars -except final_cell_id session_folder sessions case_id tt default_cd session_order
    
    load plane_all.mat
    curr_session = tt;
    plane_all_old = plane_all;
    
    n_planes = [1:1:5];
    
    clear plane_all
    
    for i=n_planes
        if i == 1
            real_roi_id{i,1} = final_cell_id{i,1}(:,curr_session);
        else
            real_roi_id{i,1} = final_cell_id{i,1}(:,curr_session) + sum(FOV_ROI(1:i-1));
        end
    end
    
    total_roi_id = vertcat(real_roi_id{:});
    
    for i=1:length(total_roi_id)
        for j=1:size(plane_all_old,2)
            plane_all{i,j} = plane_all_old{total_roi_id(i),j};
        end
        for ttt=1:4
            plane_all_decon{i,ttt} = plane_all_spk{total_roi_id(i),ttt};
        end
        for ttt=1:4
            plane_all_decon_re{i,ttt} = plane_all_spk_re{total_roi_id(i),ttt};
        end
    end
    
    mkdir(session_folder);
    cd(session_folder);
    disp(strcat('saving_',sessions{tt,1}))
    
    if length(plane_all{1,17}) == 0
        for i=1:size(plane_all_decon_re,1)
            plane_all_decon_re{i,3} = plane_all_decon_re{i,4};
            plane_all_decon{i,3} = plane_all_decon{i,4};
            plane_all{i,11} = plane_all{i,12};
            plane_all{i,13} = plane_all{i,14};
            plane_all{i,15} = plane_all{i,16};
            plane_all{i,17} = plane_all{i,18};
        end
    end
    save plane_all_matched_3sessions.mat nplanes real_roi_id total_roi_id plane_all plane_all_decon plane_all_decon_re
    
    cd(default_cd)
end

disp('done')
