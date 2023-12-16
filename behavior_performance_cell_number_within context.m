clear;clc;

currdir = dir;
default_cd = cd;
ori_cd = 'D:\Jae-Hyun Kim';

caseid = 0;
mouseid = 0;
FOVid = 0;

within_decoder_s = [];
across_decoder_s = [];
within_decoder_d = [];
across_decoder_d = [];
within_decoder_r = [];
across_decoder_r = [];

within_CD_sample = {};
across_CD_sample = {};
within_CD_delay = {};
across_CD_delay = {};
within_CD_response = {};
across_CD_response = {};

across_weight_s = {};
across_weight_d = {};
across_weight_r = {};

a1 = 0; a2 = 0;
s1 = 0; s2 = 0;
d1 = 0; d2 = 0;
r1 = 0; r2 = 0;

decoder_cri = 0;

for i=1:size(currdir,1)
    
    currfn = currdir(i).name;
    
    if length(find(currfn == '_')) == 3
        
        currfn2 = currfn;
        clear sessions
        load(currfn)
        disp(currfn)
        
        clear fncut
        fncut = find(currfn == '_');
        
        % mouse id info
        if strcmp(currfn2(6:fncut(1)-1),'35') == 1
            currmouse = strcat('BAYLORJH0',currfn2(6:fncut(1)-1));
        else
            currmouse = strcat('BAYLORJH',currfn2(6:fncut(1)-1));
        end
        
        % FOV id info
        currFOV = currfn2(fncut(1)+1:fncut(2)-1);
        
        cd(ori_cd)
        cd(currmouse)
        cd(currFOV)
        
        FOVid = FOVid + 1;
        perf_all{FOVid,1} = currmouse;
        perf_all{FOVid,2} = currFOV;
        
        % sorted out matched cell numbers
        % matched and reliable cells
        perf_all{FOVid,6}(1,1) = length(rel_cell_id);
        % matched but non-reliable cells
        perf_all{FOVid,6}(2,1) = length(non_rel_cell_id);
        % all matched cells
        perf_all{FOVid,6}(3,1) = length(rel_cell_id)+length(non_rel_cell_id);
        
        clear sessions_r
        
        for tt=1:size(sessions,1)
            sessions_r{tt,1} = convertStringsToChars(sessions{tt,1});
            sessions_r{tt,1}(find(sessions_r{tt,1} == '_')) = '-';
        end
        perf_all{FOVid,3} = sessions_r;
        
        for tt=1:size(sessions_r,1)
            for zz=1:size(sessions_r,1)
                if tt < zz
                    perf_all{FOVid,4}(tt,zz) = between(datetime(sessions_r{tt,1}),datetime(sessions_r{zz,1}));
                end
            end
        end
        
        ori_cd2 = cd;
        
        for j=1:size(sessions)
            cd(sessions{j,1})
            load plane_all.mat plane_all_spk_re
            for z=1:4
                % information about the behavior performance
                perf_all{FOVid,5}(j,z) = size(plane_all_spk_re{1,z},1);
            end
            
            % information about original number of cells
            perf_all{FOVid,7}(j,1) = size(plane_all_spk_re,1);
            
            cd(ori_cd2)
            
            clear plane_all_spk_re
        end
        
        
        cd(default_cd)
        
        elseif length(find(currfn == '_')) > 5
        
        clear cutpoint
        cutpoint = find(currfn == '_');
        
        if strcmp(currfn(cutpoint(5)+1:cutpoint(6)-1),'sample')
            load(currfn,'decoder_s_s','CD_proj_s_test','orthonormal_basis_s')
            
            clear val_session
            val_session = [];
            for t=1:size(sessions,1)
                if decoder_s_s(t,t) > decoder_cri
                    val_session = [val_session t];
                end
            end
            
            clear temp_decoder
            if length(val_session) > 1
                for t=1:length(val_session)
                    s1 = s1 + 1;
                    within_decoder_s(s1,1) = decoder_s_s(val_session(t),val_session(t));
                    
                    clear currtemp1 currtemp2
                    currtemp1 = [];
                    currtemp2 = [];
                    for zzz=1:3
                        currtemp1 = [currtemp1; CD_proj_s_test{val_session(t),val_session(t)}{zzz,2}(:,1)'];
                        currtemp2 = [currtemp2; CD_proj_s_test{val_session(t),val_session(t)}{zzz,1}(:,1)'];
                    end
                    
                    within_CD_sample{1,1}(s1,:) = mean(currtemp1,1);
                    within_CD_sample{1,2}(s1,:) = mean(currtemp2,1);
                    
                    for tt=1:length(val_session)
                        if tt > t
                            s2 = s2 + 1;
                            across_decoder_s(s2,1) = decoder_s_s(val_session(tt),val_session(t));
                            across_decoder_s(s2,2) = decoder_s_s(val_session(t),val_session(tt));
                            across_decoder_s(s2,3) = datenum(sessions{val_session(tt),1}) - datenum(sessions{val_session(t),1});
                            
                            across_weight_s{s2,1} = orthonormal_basis_s{t,1}(:,1);
                            across_weight_s{s2,2} = orthonormal_basis_s{tt,1}(:,1);
                            
                            clear currtemp1 currtemp2
                            currtemp1 = [];
                            currtemp2 = [];
                            currtemp3 = [];
                            currtemp4 = [];
                            for zzz=1:3
                                currtemp1 = [currtemp1; CD_proj_s_test{val_session(tt),val_session(t)}{zzz,2}(:,1)'];
                                currtemp2 = [currtemp2; CD_proj_s_test{val_session(tt),val_session(t)}{zzz,1}(:,1)'];
                                currtemp3 = [currtemp3; CD_proj_s_test{val_session(t),val_session(tt)}{zzz,2}(:,1)'];
                                currtemp4 = [currtemp4; CD_proj_s_test{val_session(t),val_session(tt)}{zzz,1}(:,1)'];
                            end
                            
                            across_CD_sample{1,1}(s2,:) = mean(currtemp1,1);
                            across_CD_sample{1,2}(s2,:) = mean(currtemp2,1);
                            across_CD_sample{1,3}(s2,:) = mean(currtemp3,1);
                            across_CD_sample{1,4}(s2,:) = mean(currtemp4,1);
                        end
                    end
                end
            end
            
        elseif strcmp(currfn(cutpoint(5)+1:cutpoint(6)-1),'delay')
            
            load(currfn,'decoder_d_d','CD_proj_d_test','orthonormal_basis_d')
            
            clear val_session
            val_session = [];
            for t=1:size(sessions,1)
                if decoder_d_d(t,t) > decoder_cri
                    val_session = [val_session t];
                end
            end
            
            clear temp_decoder
            if length(val_session) > 1
                for t=1:length(val_session)
                    d1 = d1 + 1;
                    within_decoder_d(d1,1) = decoder_d_d(val_session(t),val_session(t));
                    
                    clear currtemp1 currtemp2
                    currtemp1 = [];
                    currtemp2 = [];
                    for zzz=1:3
                        currtemp1 = [currtemp1; CD_proj_d_test{val_session(t),val_session(t)}{zzz,2}(:,1)'];
                        currtemp2 = [currtemp2; CD_proj_d_test{val_session(t),val_session(t)}{zzz,1}(:,1)'];
                    end
                    
                    within_CD_delay{1,1}(d1,:) = mean(currtemp1,1);
                    within_CD_delay{1,2}(d1,:) = mean(currtemp2,1);
                    
                    for tt=1:length(val_session)
                        if tt > t
                            d2 = d2 + 1;
                            across_decoder_d(d2,1) = decoder_d_d(val_session(tt),val_session(t));
                            across_decoder_d(d2,2) = decoder_d_d(val_session(t),val_session(tt));
                            across_decoder_d(d2,3) = datenum(sessions{val_session(tt),1}) - datenum(sessions{val_session(t),1});
                            
                            across_weight_d{d2,1} = orthonormal_basis_d{t,1}(:,1);
                            across_weight_d{d2,2} = orthonormal_basis_d{tt,1}(:,1);
                            
                            clear currtemp1 currtemp2
                            currtemp1 = [];
                            currtemp2 = [];
                            currtemp3 = [];
                            currtemp4 = [];
                            for zzz=1:3
                                currtemp1 = [currtemp1; CD_proj_d_test{val_session(tt),val_session(t)}{zzz,2}(:,1)'];
                                currtemp2 = [currtemp2; CD_proj_d_test{val_session(tt),val_session(t)}{zzz,1}(:,1)'];
                                currtemp3 = [currtemp3; CD_proj_d_test{val_session(t),val_session(tt)}{zzz,2}(:,1)'];
                                currtemp4 = [currtemp4; CD_proj_d_test{val_session(t),val_session(tt)}{zzz,1}(:,1)'];
                            end
                            
                            across_CD_delay{1,1}(d2,:) = mean(currtemp1,1);
                            across_CD_delay{1,2}(d2,:) = mean(currtemp2,1);
                            across_CD_delay{1,3}(d2,:) = mean(currtemp3,1);
                            across_CD_delay{1,4}(d2,:) = mean(currtemp4,1);
                        end
                    end
                end
            end
            
        elseif strcmp(currfn(cutpoint(5)+1:cutpoint(6)-1),'response')
            
            load(currfn,'decoder_r_r','CD_proj_r_test','orthonormal_basis_r')
            
            clear val_session
            val_session = [];
            for t=1:size(sessions,1)
                if decoder_r_r(t,t) > decoder_cri
                    val_session = [val_session t];
                end
            end
            
            clear temp_decoder
            if length(val_session) > 1
                for t=1:length(val_session)
                    r1 = r1 + 1;
                    within_decoder_r(r1,1) = decoder_r_r(val_session(t),val_session(t));
                    
                    clear currtemp1 currtemp2
                    currtemp1 = [];
                    currtemp2 = [];
                    for zzz=1:3
                        currtemp1 = [currtemp1; CD_proj_r_test{val_session(t),val_session(t)}{zzz,2}(:,1)'];
                        currtemp2 = [currtemp2; CD_proj_r_test{val_session(t),val_session(t)}{zzz,1}(:,1)'];
                    end
                    
                    within_CD_response{1,1}(r1,:) = mean(currtemp1,1);
                    within_CD_response{1,2}(r1,:) = mean(currtemp2,1);
                    
                    for tt=1:length(val_session)
                        if tt > t
                            r2 = r2 + 1;
                            across_decoder_r(r2,1) = decoder_r_r(val_session(tt),val_session(t));
                            across_decoder_r(r2,2) = decoder_r_r(val_session(t),val_session(tt));
                            across_decoder_r(r2,3) = datenum(sessions{val_session(tt),1}) - datenum(sessions{val_session(t),1});
                            
                            across_weight_r{r2,1} = orthonormal_basis_r{t,1}(:,1);
                            across_weight_r{r2,2} = orthonormal_basis_r{tt,1}(:,1);
                            
                            clear currtemp1 currtemp2
                            currtemp1 = [];
                            currtemp2 = [];
                            currtemp3 = [];
                            currtemp4 = [];
                            for zzz=1:3
                                currtemp1 = [currtemp1; CD_proj_r_test{val_session(tt),val_session(t)}{zzz,2}(:,1)'];
                                currtemp2 = [currtemp2; CD_proj_r_test{val_session(tt),val_session(t)}{zzz,1}(:,1)'];
                                currtemp3 = [currtemp3; CD_proj_r_test{val_session(t),val_session(tt)}{zzz,2}(:,1)'];
                                currtemp4 = [currtemp4; CD_proj_r_test{val_session(t),val_session(tt)}{zzz,1}(:,1)'];
                            end
                            
                            across_CD_response{1,1}(r2,:) = mean(currtemp1,1);
                            across_CD_response{1,2}(r2,:) = mean(currtemp2,1);
                            across_CD_response{1,3}(r2,:) = mean(currtemp1,1);
                            across_CD_response{1,4}(r2,:) = mean(currtemp2,1);
                        end
                    end
                end
            end
            
            
        end
    end
end


save behavior_cell_number.mat perf_all

%%

clear;clc; 
load behavior_cell_number.mat

mouseid = 1;
final_all = {};
for i=1:size(perf_all,1)
    
    currmouseid = perf_all{i,1};
    if i == 1
        premouseid = perf_all{i,1};
        final_all{mouseid,1} = perf_all{i,3};
        final_all{mouseid,2} = perf_all{i,5};
        
        final_all{mouseid,3} = currmouseid;
    else
        premouseid = perf_all{i-1,1};
        if strcmp(currmouseid,premouseid) ~= 1
            mouseid = mouseid + 1;
            final_all{mouseid,1} = perf_all{i,3};
            final_all{mouseid,2} = perf_all{i,5};
            final_all{mouseid,3} = currmouseid;
        else
            
            final_all{mouseid,1} = vertcat(final_all{mouseid,1},perf_all{i,3});
            final_all{mouseid,2} = vertcat(final_all{mouseid,2},perf_all{i,5});
            
        end
    end
    
end


% plot the behavior performance
sz = 40;
sz2 = 15;
limcri = 61;
figure
%colorcode = rand(size(final_all,1),3);
colorcode = hsv(size(final_all,1));

for i=1:size(final_all,1)
    
    for j=1:size(final_all{i,2},1)
        final_all{i,5}(j,1) = sum(final_all{i,2}(j,1:2))/sum(final_all{i,2}(j,:));
    end
    
    clear aa index
    [aa index] = sortrows(final_all{i,1});
    final_all{i,4} = final_all{i,1}(index);
    for j=1:length(final_all{i,4})
        final_all{i,7}(j,1) = datenum(final_all{i,4}(j,1)) - datenum(final_all{i,4}(1,1));
    end
    final_all{i,6} = final_all{i,5}(index);
    
    %
    final_all_re{i,7} = final_all{i,7};
    final_all_re{i,7}(find(final_all{i,7} > limcri)) = limcri-1;
    
    hold on
    scatter(final_all_re{i,7},final_all{i,6},sz,colorcode(i,:),'filled')
    hold on
    for j=1:length(final_all{i,7})-1
        line([final_all_re{i,7}(j) final_all_re{i,7}(j+1)],[final_all{i,6}(j) final_all{i,6}(j+1)],'color',colorcode(i,:))     
    end

end

session_no = 0;
perf_comb = [];
for j=1:size(final_all,1)
    session_no = session_no + length(final_all{j,5});
    perf_comb = vertcat(perf_comb,final_all{j,5});
end

ylim([0 1])
xlim([-5 limcri])
xlabel('Delta days','fontsize',15)
ylabel('Fraction correct','fontsize',15)
hold on
title(strcat(num2str(size(final_all,1)),' mice, ',' ',num2str(session_no),' sessions'),'fontsize',20)

ax=gca;
ax.XAxis.FontSize=13;
ax.YAxis.FontSize=13;
   
set(gcf,'color','w')

default_cd = pwd;
cd('figures')
exportgraphics(gcf,'withincontextovertime_behavior_performance.emf','ContentType','vector')
cd(default_cd)

clearvars -except colorcode
%%
clc; 
load behavior_cell_number.mat

% extracting mouse id (number)
for i=1:size(perf_all,1)
    mouseid(i,1) = str2num(perf_all{i,1}(9:end));
end

clear final_all
mousepool = unique(mouseid);
for i=1:length(mousepool)
    
    clear currmouseFOVs
    currmouseFOVs = find(mouseid == mousepool(i));
    
    % mouse id
    final_all{i,1} = perf_all{currmouseFOVs(1),1};
    
    % performance and delta days
    for j=1:length(currmouseFOVs)
        clear temp
        temp = perf_all{currmouseFOVs(j),5};
        % early session performance
        final_all{i,2}(j,1) = sum(temp(1,1:2))/sum(temp(1,1:4));
        % late session performance
        final_all{i,3}(j,1) = sum(temp(2,1:2))/sum(temp(2,1:4));
        
        % delta days
        final_all{i,4}(j,1) = datenum(perf_all{currmouseFOVs(j),3}(2)) - datenum(perf_all{currmouseFOVs(j),3}(1));
        
        % number of cells from early session
        final_all{i,5}(j,1) = perf_all{currmouseFOVs(j),7}(1);
        
        % number of cells from late session
        final_all{i,6}(j,1) = perf_all{currmouseFOVs(j),7}(2);
        
        % number of all matched cells
        final_all{i,7}(j,1) = perf_all{currmouseFOVs(j),6}(3);
        
        % number of matched and reliable cells
        final_all{i,8}(j,1) = perf_all{currmouseFOVs(j),6}(1);
        
        % number of matched and non-reliable cells
        final_all{i,9}(j,1) = perf_all{currmouseFOVs(j),6}(2);
    end
end


% original number of cells, matched number of cells, reliable/non-reliable matched cells

sz = 20;
figure

cell_no = [];
for i=1:size(final_all,1)
    clear temp
    for j=1:5
        temp(:,j) = final_all{i,j+4};
    end
    
    cell_no = vertcat(cell_no,temp);
end

for i=1:2
    subplot(1,2,i)
    hold on
    bar(mean(cell_no,1),.6,'FaceColor','k','FaceAlpha',.1)
    %errorbar([1:1:5],mean(cell_no,1),std(cell_no,1)/sqrt(size(cell_no,1)),'k','LineWidth',1)
    errorbar([1:1:5],mean(cell_no,1),std(cell_no,1),'k','LineWidth',1)
end

for i=1:size(final_all,1)
    
    % plot with lines
    subplot(1,2,1)
    for j=1:5
        % early session cell number
        hold on
        scatter(j*ones(1,length(final_all{i,4+j})),final_all{i,4+j},sz,colorcode(i,:),'filled','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
    end
    
    hold on
    for j=1:length(final_all{i,5})
        for k=1:4
            line([0+k 1+k],[final_all{i,4+k}(j) final_all{i,5+k}(j)],'color',[colorcode(i,:) .3],'LineWidth',.1)
        end
    end
    
    % plot without lines
    subplot(1,2,2)
    for j=1:5
        % early session cell number
        hold on
        scatter(j*ones(1,length(final_all{i,4+j})),final_all{i,4+j},sz,colorcode(i,:),'filled','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
    end
    
end

for i=1:2
    subplot(1,2,i)
    ylim([0 4500])
    xlim([0 6])
    xticks([1:1:5])
    xticklabels({'Context 1','Context 2','Matched','Matched and rel.','Matched and non-rel.'})
    xtickangle(45)
    ylabel('Cell number','fontsize',15)
end

set(gcf,'color','w')
%%
default_cd = pwd;
cd('figures')
exportgraphics(gcf,'withincontextovertime_cell_numbers.emf','ContentType','vector')
cd(default_cd)








