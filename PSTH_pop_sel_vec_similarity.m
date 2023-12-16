clear; clc;

currdir = dir;

stat_p_all = [];
stat_sel_all = [];
resp_all_b = [];
resp_all_s = [];
resp_all_d = [];
resp_all_r = [];

within_psth_sim = [];
across_psth_sim = [];
within_popsel_sim = [];
across_popsel_sim = [];
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
s3 = 0;
d3 = 0;
r3 = 0;
s_accum = [];
d_accum = [];
r_accum = [];

decoder_cri = .695;

for i=1:size(currdir,1)
    
    currfn = currdir(i).name;
    
    if length(find(currfn == '_')) == 3
        
        currfn2 = currfn;
        
        load(currfn)
        disp(currfn)
        
        if size(sessions,1) == 3
            stat_p = stat_p(:,[1 2 4 5 7 8]);
            stat_sel = stat_sel(:,[1 2 4 5 7 8]);            
        elseif size(sessions,1) == 4
            stat_p = stat_p(:,[1 2 5 6 9 10]);
            stat_sel = stat_sel(:,[1 2 5 6 9 10]);
        elseif size(sessions,1) == 5
            stat_p = stat_p(:,[1 2 6 7 11 12]);
            stat_sel = stat_sel(:,[1 2 6 7 11 12]);
        end
        
        stat_p_all = vertcat(stat_p_all,stat_p(rel_cell_id,:));
        stat_sel_all = vertcat(stat_sel_all,stat_sel(rel_cell_id,:));
            
        resp_all_b = vertcat(resp_all_b,horzcat((resp_all{1,1}(rel_cell_id,2)-resp_all{1,1}(rel_cell_id,1))./(resp_all{1,1}(rel_cell_id,2)+resp_all{1,1}(rel_cell_id,1)),...
            (resp_all{2,1}(rel_cell_id,2)-resp_all{2,1}(rel_cell_id,1))./(resp_all{2,1}(rel_cell_id,2)+resp_all{2,1}(rel_cell_id,1))));
        resp_all_s = vertcat(resp_all_s,horzcat((resp_all{1,2}(rel_cell_id,2)-resp_all{1,2}(rel_cell_id,1))./(resp_all{1,2}(rel_cell_id,2)+resp_all{1,2}(rel_cell_id,1)),...
            (resp_all{2,2}(rel_cell_id,2)-resp_all{2,2}(rel_cell_id,1))./(resp_all{2,2}(rel_cell_id,2)+resp_all{2,2}(rel_cell_id,1))));
        resp_all_d = vertcat(resp_all_d,horzcat((resp_all{1,3}(rel_cell_id,2)-resp_all{1,3}(rel_cell_id,1))./(resp_all{1,3}(rel_cell_id,2)+resp_all{1,3}(rel_cell_id,1)),...
            (resp_all{2,3}(rel_cell_id,2)-resp_all{2,3}(rel_cell_id,1))./(resp_all{2,3}(rel_cell_id,2)+resp_all{2,3}(rel_cell_id,1))));
        resp_all_r = vertcat(resp_all_r,horzcat((resp_all{1,4}(rel_cell_id,2)-resp_all{1,4}(rel_cell_id,1))./(resp_all{1,4}(rel_cell_id,2)+resp_all{1,4}(rel_cell_id,1)),...
            (resp_all{2,4}(rel_cell_id,2)-resp_all{2,4}(rel_cell_id,1))./(resp_all{2,4}(rel_cell_id,2)+resp_all{2,4}(rel_cell_id,1))));
        
        clear cutpoint
        cutpoint = find(currfn == '_');
        mouseid = str2num(currfn(6:cutpoint(1)-1));
        
        for t=1:size(sessions,1)
            a1 = a1 + 1;
            within_psth_sim(a1,1) = take_home{1,1}(1,t);
            within_psth_sim(a1,2) = take_home{1,2}(1,t);
            within_psth_sim(a1,3) = mouseid;
            within_popsel_sim(a1,1) = take_home{2,1}(1,t);
            within_popsel_sim(a1,2) = take_home{2,2}(1,t);
            within_popsel_sim(a1,3) = mouseid;
            
            for tt=1:size(sessions,1)
                if tt > t
                    a2 = a2 + 1;
                    across_psth_sim(a2,1) = take_home{1,1}(tt,t);
                    across_psth_sim(a2,2) = take_home{1,2}(tt,t);
                    across_psth_sim(a2,3) = datenum(sessions{tt,1}) - datenum(sessions{t,1});
                    across_psth_sim(a2,4) = mouseid;
                    
                    across_popsel_sim(a2,1) = take_home{2,1}(tt,t);
                    across_popsel_sim(a2,2) = take_home{2,2}(tt,t);
                    across_popsel_sim(a2,3) = datenum(sessions{tt,1}) - datenum(sessions{t,1});
                    across_popsel_sim(a2,4) = mouseid;
                end
            end
        end
        
    elseif length(find(currfn == '_')) > 5
        
        clear cutpoint
        cutpoint = find(currfn == '_');
        
        if strcmp(currfn(cutpoint(5)+1:cutpoint(6)-1),'sample')
            load(currfn,'decoder_s_s','CD_proj_s_test','orthonormal_basis_s')
            
            mouseid = str2num(currfn(6:cutpoint(1)-1));
            FOV = str2num(currfn(cutpoint(1)+4:cutpoint(2)-1));
            
            clear val_session
            val_session = [];
            for t=1:size(sessions,1)
                if decoder_s_s(t,t) > decoder_cri
                    val_session = [val_session t];
                end
            end
            
            clear temp_decoder
            if length(val_session) > 1
                s3 = s3 + 1;
                s_accum(s3,1) = mouseid;
                s_accum(s3,2) = FOV;
                s_accum(s3,3) = length(val_session);
                
                for t=1:length(val_session)
                    s1 = s1 + 1;
                    within_decoder_s(s1,1) = decoder_s_s(val_session(t),val_session(t));
                    within_decoder_s(s1,2) = mouseid;
                    within_decoder_s(s1,3) = FOV;
                    
                    clear currtemp1 currtemp2
                    currtemp1 = [];
                    currtemp2 = [];
                    for zzz=1:2
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
                            across_decoder_s(s2,4) = mouseid;
                            across_decoder_s(s2,5) = decoder_s_s(val_session(t),val_session(t));
                            across_decoder_s(s2,6) = decoder_s_s(val_session(tt),val_session(tt));
                            
                            across_weight_s{s2,1} = orthonormal_basis_s{t,1}(:,1);
                            across_weight_s{s2,2} = orthonormal_basis_s{tt,1}(:,1);
                            
                            clear currtemp1 currtemp2
                            currtemp1 = [];
                            currtemp2 = [];
                            currtemp3 = [];
                            currtemp4 = [];
                            for zzz=1:2
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
            mouseid = str2num(currfn(6:cutpoint(1)-1));
            FOV = str2num(currfn(cutpoint(1)+4:cutpoint(2)-1));
            
            clear val_session
            val_session = [];
            for t=1:size(sessions,1)
                if decoder_d_d(t,t) > decoder_cri
                    val_session = [val_session t];
                end
            end
            
            clear temp_decoder
            if length(val_session) > 1
                d3 = d3 + 1;
                d_accum(d3,1) = mouseid;
                d_accum(d3,2) = FOV;
                d_accum(d3,3) = length(val_session);
                
                for t=1:length(val_session)
                    d1 = d1 + 1;
                    within_decoder_d(d1,1) = decoder_d_d(val_session(t),val_session(t));
                    within_decoder_d(d1,2) = mouseid;
                    within_decoder_d(d1,3) = FOV;
                    
                    clear currtemp1 currtemp2
                    currtemp1 = [];
                    currtemp2 = [];
                    for zzz=1:2
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
                            across_decoder_d(d2,4) = mouseid;
                            across_decoder_d(d2,5) = decoder_d_d(val_session(t),val_session(t));
                            across_decoder_d(d2,6) = decoder_d_d(val_session(tt),val_session(tt));
                            
                            across_weight_d{d2,1} = orthonormal_basis_d{t,1}(:,1);
                            across_weight_d{d2,2} = orthonormal_basis_d{tt,1}(:,1);
                            
                            clear currtemp1 currtemp2
                            currtemp1 = [];
                            currtemp2 = [];
                            currtemp3 = [];
                            currtemp4 = [];
                            for zzz=1:2
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
            mouseid = str2num(currfn(6:cutpoint(1)-1));
            FOV = str2num(currfn(cutpoint(1)+4:cutpoint(2)-1));
            
            clear val_session
            val_session = [];
            for t=1:size(sessions,1)
                if decoder_r_r(t,t) > decoder_cri
                    val_session = [val_session t];
                end
            end
            
            clear temp_decoder
            if length(val_session) > 1
                r3 = r3 + 1;
                r_accum(r3,1) = mouseid;
                r_accum(r3,2) = FOV;
                r_accum(r3,3) = length(val_session);
                
                for t=1:length(val_session)
                    r1 = r1 + 1;
                    within_decoder_r(r1,1) = decoder_r_r(val_session(t),val_session(t));
                    within_decoder_r(r1,2) = mouseid;
                    within_decoder_r(r1,3) = FOV;
                     
                    clear currtemp1 currtemp2
                    currtemp1 = [];
                    currtemp2 = [];
                    for zzz=1:2
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
                            across_decoder_r(r2,4) = mouseid;
                            across_decoder_r(r2,5) = decoder_r_r(val_session(t),val_session(t));
                            across_decoder_r(r2,6) = decoder_r_r(val_session(tt),val_session(tt));
                            
                            across_weight_r{r2,1} = orthonormal_basis_r{t,1}(:,1);
                            across_weight_r{r2,2} = orthonormal_basis_r{tt,1}(:,1);
                            
                            clear currtemp1 currtemp2
                            currtemp1 = [];
                            currtemp2 = [];
                            currtemp3 = [];
                            currtemp4 = [];
                            for zzz=1:2
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
%%

decoder_s_re = {};
decoder_d_re = {};
decoder_r_re = {};
mouseid_pool = {'117','119','122','130','133','134','160','35'};
for i=1:length(mouseid_pool)
%     clear temp
%     temp = find(within_decoder_s(:,2) == str2num(cell2mat(mouseid_pool(i))));
%     decoder_s_re{i,1} = within_decoder_s(temp,1);
    
    clear temp
    temp = find(across_decoder_s(:,4) == str2num(cell2mat(mouseid_pool(i))));
    decoder_s_re{i,1} = across_decoder_s(temp,3);   % delta days
    decoder_s_re{i,2} = across_decoder_s(temp,5);   % train early test early
    decoder_s_re{i,3} = across_decoder_s(temp,6);   % train late test late
    decoder_s_re{i,4} = across_decoder_s(temp,2);   % train early test late
    decoder_s_re{i,5} = across_decoder_s(temp,1);   % train late test early
    decoder_s_re{i,6} = mean(horzcat(decoder_s_re{i,2},decoder_s_re{i,3}),2);   % within session
    decoder_s_re{i,7} = mean(horzcat(decoder_s_re{i,4},decoder_s_re{i,5}),2);   % across sessions
    
%     clear temp
%     temp = find(within_decoder_d(:,2) == str2num(cell2mat(mouseid_pool(i))));
%     decoder_d_re{i,1} = within_decoder_d(temp,1);
    
    clear temp
    temp = find(across_decoder_d(:,4) == str2num(cell2mat(mouseid_pool(i))));
    decoder_d_re{i,1} = across_decoder_d(temp,3);   % delta days
    decoder_d_re{i,2} = across_decoder_d(temp,5);   % train early test early
    decoder_d_re{i,3} = across_decoder_d(temp,6);   % train late test late
    decoder_d_re{i,4} = across_decoder_d(temp,2);   % train early test late
    decoder_d_re{i,5} = across_decoder_d(temp,1);   % train late test early
    decoder_d_re{i,6} = mean(horzcat(decoder_d_re{i,2},decoder_d_re{i,3}),2);   % within session
    decoder_d_re{i,7} = mean(horzcat(decoder_d_re{i,4},decoder_d_re{i,5}),2);   % across sessions
    
    
%     clear temp
%     temp = find(within_decoder_r(:,2) == str2num(cell2mat(mouseid_pool(i))));
%     decoder_r_re{i,1} = within_decoder_r(temp,1);
    
    clear temp
    temp = find(across_decoder_r(:,4) == str2num(cell2mat(mouseid_pool(i))));
    decoder_r_re{i,1} = across_decoder_r(temp,3);   % delta days
    decoder_r_re{i,2} = across_decoder_r(temp,5);   % train early test early
    decoder_r_re{i,3} = across_decoder_r(temp,6);   % train late test late
    decoder_r_re{i,4} = across_decoder_r(temp,2);   % train early test late
    decoder_r_re{i,5} = across_decoder_r(temp,1);   % train late test early
    decoder_r_re{i,6} = mean(horzcat(decoder_r_re{i,2},decoder_r_re{i,3}),2);   % within session
    decoder_r_re{i,7} = mean(horzcat(decoder_r_re{i,4},decoder_r_re{i,5}),2);   % across sessions
end
%%
sz = 25;
colorcode = hsv(size(decoder_s_re,1));
figure
for i=1:size(colorcode,1)
    subplot(3,2,1)
    hold on
    scatter(rand(size(decoder_s_re{i,2},1),1)*-5,decoder_s_re{i,2},sz,colorcode(i,:),'filled')
    hold on
    scatter(rand(size(decoder_s_re{i,3},1),1)*-5,decoder_s_re{i,3},sz,colorcode(i,:),'filled')
    hold on
    scatter(decoder_s_re{i,1},decoder_s_re{i,4},sz,colorcode(i,:),'filled')
    hold on
    scatter(decoder_s_re{i,1},decoder_s_re{i,5},sz,colorcode(i,:),'filled')

    subplot(3,2,3)
    hold on
    scatter(rand(size(decoder_d_re{i,2},1),1)*-5,decoder_d_re{i,2},sz,colorcode(i,:),'filled')
    hold on
    scatter(rand(size(decoder_d_re{i,3},1),1)*-5,decoder_d_re{i,3},sz,colorcode(i,:),'filled')
    hold on
    scatter(decoder_d_re{i,1},decoder_d_re{i,4},sz,colorcode(i,:),'filled')
    hold on
    scatter(decoder_d_re{i,1},decoder_d_re{i,5},sz,colorcode(i,:),'filled')
    
    subplot(3,2,5)
    hold on
    scatter(rand(size(decoder_r_re{i,2},1),1)*-5,decoder_r_re{i,2},sz,colorcode(i,:),'filled')
    hold on
    scatter(rand(size(decoder_r_re{i,3},1),1)*-5,decoder_r_re{i,3},sz,colorcode(i,:),'filled')
    hold on
    scatter(decoder_r_re{i,1},decoder_r_re{i,4},sz,colorcode(i,:),'filled')
    hold on
    scatter(decoder_r_re{i,1},decoder_r_re{i,5},sz,colorcode(i,:),'filled')
    
    
    subplot(3,2,2)
    hold on
    scatter(rand(size(decoder_s_re{i,6},1),1)*-5,decoder_s_re{i,6},sz,colorcode(i,:),'filled')
    hold on
    scatter(decoder_s_re{i,1},decoder_s_re{i,7},sz,colorcode(i,:),'filled')

    subplot(3,2,4)
    hold on
    scatter(rand(size(decoder_d_re{i,6},1),1)*-5,decoder_d_re{i,6},sz,colorcode(i,:),'filled')
    hold on
    scatter(decoder_d_re{i,1},decoder_d_re{i,7},sz,colorcode(i,:),'filled')
    
    subplot(3,2,6)
    hold on
    scatter(rand(size(decoder_r_re{i,6},1),1)*-5,decoder_r_re{i,6},sz,colorcode(i,:),'filled')
    hold on
    scatter(decoder_r_re{i,1},decoder_r_re{i,7},sz,colorcode(i,:),'filled')
end

for i=1:6
    subplot(3,2,i)
    ylim([0 1])
    ylabel('Decoding accuracy','fontsize',15)
    xlim([-8 62])
    xlabel('Delta days','fontsize',15)
    hold on
    line([0 62], [.5 .5],'color','k','LineStyle',':')
    ax=gca;
    ax.XAxis.FontSize=13;
    ax.YAxis.FontSize=13;
end


subplot(3,2,1)
hold on
title(strcat('8 mice,',num2str(round(length(within_decoder_s)/2)),' FOVs'),'fontsize',15)

subplot(3,2,3)
hold on
title(strcat('8 mice,',num2str(round(length(within_decoder_d)/2)),' FOVs'),'fontsize',15)

subplot(3,2,5)
hold on
title(strcat('8 mice,',num2str(round(length(within_decoder_r)/2)),' FOVs'),'fontsize',15)



subplot(3,2,2)
hold on
title(strcat('8 mice, average, ',num2str(round(length(within_decoder_s)/2)),' FOVs'),'fontsize',15)

subplot(3,2,4)
hold on
title(strcat('8 mice, avearge, ',num2str(round(length(within_decoder_d)/2)),' FOVs'),'fontsize',15)

subplot(3,2,6)
hold on
title(strcat('8 mice, average, ',num2str(round(length(within_decoder_r)/2)),' FOVs'),'fontsize',15)

set(gcf,'color','w')

%exportgraphics(gcf,'SameContextDAovertime_SDR.emf','ContentType','vector')
%%






sz = 7;
figure
subplot(3,3,1)
hold on
scatter(rand(1,s1)*-5,within_decoder_s(:,1),sz,'k')
hold on
scatter(across_decoder_s(:,3),across_decoder_s(:,1),sz,'b')
ylim([0 1])
ylabel('accuracy')
xlim([-8 48])
xlabel('delta days')
title('sample: late -> early')

subplot(3,3,2)
hold on
scatter(rand(1,s1)*-5,within_decoder_s(:,1),sz,'k')
hold on
scatter(across_decoder_s(:,3),across_decoder_s(:,2),sz,'b')
ylim([0 1])
ylabel('accuracy')
xlim([-8 48])
xlabel('delta days')
title('sample: early -> late')

subplot(3,3,3)
hold on
scatter(across_decoder_s(:,2),across_decoder_s(:,1),sz,'k')
ylim([.4 1])
xlim([.4 1])
ylabel('late -> early')
xlabel('early -> late')
hold on
title(strcat('p=',num2str(signrank(across_decoder_s(:,2),across_decoder_s(:,1)))))

subplot(3,3,4)
hold on
scatter(rand(1,d1)*-5,within_decoder_d(:,1),sz,'k')
hold on
scatter(across_decoder_d(:,3),across_decoder_d(:,1),sz,'b')
ylim([0 1])
ylabel('accuracy')
xlim([-8 62])
xlabel('delta days')
title('delay: late -> early')

subplot(3,3,5)
hold on
scatter(rand(1,d1)*-5,within_decoder_d(:,1),sz,'k')
hold on
scatter(across_decoder_d(:,3),across_decoder_d(:,2),sz,'b')
ylim([0 1])
ylabel('accuracy')
xlim([-8 62])
xlabel('delta days')
title('delay: early -> late')

subplot(3,3,6)
hold on
scatter(across_decoder_d(:,2),across_decoder_d(:,1),sz,'k')
ylim([.4 1])
xlim([.4 1])
ylabel('late -> early')
xlabel('early -> late')
hold on
title(strcat('p=',num2str(signrank(across_decoder_d(:,2),across_decoder_d(:,1)))))

subplot(3,3,7)
hold on
scatter(rand(1,r1)*-5,within_decoder_r(:,1),sz,'k')
hold on
scatter(across_decoder_r(:,3),across_decoder_r(:,1),sz,'b')
ylim([0 1])
ylabel('accuracy')
xlim([-8 62])
xlabel('delta days')
title('response: late -> early')

subplot(3,3,8)
hold on
scatter(rand(1,r1)*-5,within_decoder_r(:,1),sz,'k')
hold on
scatter(across_decoder_r(:,3),across_decoder_r(:,2),sz,'b')
ylim([0 1])
ylabel('accuracy')
xlim([-8 62])
xlabel('delta days')
hold on
title('response: early -> late')

subplot(3,3,9)
hold on
scatter(across_decoder_r(:,2),across_decoder_r(:,1),sz,'k')
ylim([.4 1])
xlim([.4 1])
ylabel('late -> early')
xlabel('early -> late')
hold on
title(strcat('p=',num2str(signrank(across_decoder_r(:,2),across_decoder_r(:,1)))))

hold on
sgtitle(strcat('within session decoder accuracy >',num2str(decoder_cri)))

binloc = [1:1:47];
bintick = round([1.57 2.87 4.17]*6);
xlimbd = round([1.07 6.17]*6);
xlabeltext = [-2.6, -1.3, 0];

%%
figure
subplot(3,3,1)
hold on
plot(within_CD_sample{1,1}','color',[1 0 0 .1],'LineWidth',.1)
hold on
plot(within_CD_sample{1,2}','color',[0 0 1 .1],'LineWidth',.1)
hold on
plot(mean(within_CD_sample{1,1},1),'color',[1 0 0],'LineWidth',2)
hold on
plot(mean(within_CD_sample{1,2},1),'color',[0 0 1],'LineWidth',2)
ylim([-4 6])

xlim(xlimbd)
xticks(bintick)
xticklabels(xlabeltext)
for tt=1:length(bintick)
    line([bintick(tt) bintick(tt)], [-2 4], 'color','k','LineStyle',':')
end
xlabel('time from go cue(s)')
ylabel('CD projection(a.u.)')

subplot(3,3,2)
hold on
plot(across_CD_sample{1,1}','color',[1 0 0 .1],'LineWidth',.1)
hold on
plot(across_CD_sample{1,2}','color',[0 0 1 .1],'LineWidth',.1)
hold on
plot(mean(across_CD_sample{1,1},1),'color',[1 0 0],'LineWidth',2)
hold on
plot(mean(across_CD_sample{1,2},1),'color',[0 0 1],'LineWidth',2)
ylim([-4 6])

xlim(xlimbd)
xticks(bintick)
xticklabels(xlabeltext)
for tt=1:length(bintick)
    line([bintick(tt) bintick(tt)], [-2 4], 'color','k','LineStyle',':')
end
xlabel('time from go cue(s)')
ylabel('CD projection(a.u.)')

subplot(3,3,3)
hold on
plot(across_CD_sample{1,3}','color',[1 0 0 .1],'LineWidth',.1)
hold on
plot(across_CD_sample{1,4}','color',[0 0 1 .1],'LineWidth',.1)
hold on
plot(mean(across_CD_sample{1,3},1),'color',[1 0 0],'LineWidth',2)
hold on
plot(mean(across_CD_sample{1,4},1),'color',[0 0 1],'LineWidth',2)
ylim([-4 6])

xlim(xlimbd)
xticks(bintick)
xticklabels(xlabeltext)
for tt=1:length(bintick)
    line([bintick(tt) bintick(tt)], [-2 4], 'color','k','LineStyle',':')
end
xlabel('time from go cue(s)')
ylabel('CD projection(a.u.)')

subplot(3,3,4)
hold on
plot(within_CD_delay{1,1}','color',[1 0 0 .1],'LineWidth',.1)
hold on
plot(within_CD_delay{1,2}','color',[0 0 1 .1],'LineWidth',.1)
hold on
plot(mean(within_CD_delay{1,1},1),'color',[1 0 0],'LineWidth',2)
hold on
plot(mean(within_CD_delay{1,2},1),'color',[0 0 1],'LineWidth',2)
ylim([-2 4])

xlim(xlimbd)
xticks(bintick)
xticklabels(xlabeltext)
for tt=1:length(bintick)
    line([bintick(tt) bintick(tt)], [-2 4], 'color','k','LineStyle',':')
end
xlabel('time from go cue(s)')
ylabel('CD projection(a.u.)')

subplot(3,3,5)
hold on
plot(across_CD_delay{1,1}','color',[1 0 0 .1],'LineWidth',.1)
hold on
plot(across_CD_delay{1,2}','color',[0 0 1 .1],'LineWidth',.1)
hold on
plot(mean(across_CD_delay{1,1},1),'color',[1 0 0],'LineWidth',2)
hold on
plot(mean(across_CD_delay{1,2},1),'color',[0 0 1],'LineWidth',2)
ylim([-2 4])

xlim(xlimbd)
xticks(bintick)
xticklabels(xlabeltext)
for tt=1:length(bintick)
    line([bintick(tt) bintick(tt)], [-2 4], 'color','k','LineStyle',':')
end
xlabel('time from go cue(s)')
ylabel('CD projection(a.u.)')

subplot(3,3,6)
hold on
plot(across_CD_delay{1,3}','color',[1 0 0 .1],'LineWidth',.1)
hold on
plot(across_CD_delay{1,4}','color',[0 0 1 .1],'LineWidth',.1)
hold on
plot(mean(across_CD_delay{1,3},1),'color',[1 0 0],'LineWidth',2)
hold on
plot(mean(across_CD_delay{1,4},1),'color',[0 0 1],'LineWidth',2)
ylim([-2 4])

xlim(xlimbd)
xticks(bintick)
xticklabels(xlabeltext)
for tt=1:length(bintick)
    line([bintick(tt) bintick(tt)], [-2 4], 'color','k','LineStyle',':')
end
xlabel('time from go cue(s)')
ylabel('CD projection(a.u.)')

subplot(3,3,7)
hold on
plot(within_CD_response{1,1}','color',[1 0 0 .1],'LineWidth',.1)
hold on
plot(within_CD_response{1,2}','color',[0 0 1 .1],'LineWidth',.1)
hold on
plot(mean(within_CD_response{1,1},1),'color',[1 0 0],'LineWidth',2)
hold on
plot(mean(within_CD_response{1,2},1),'color',[0 0 1],'LineWidth',2)
ylim([-2 4])

xlim(xlimbd)
xticks(bintick)
xticklabels(xlabeltext)
for tt=1:length(bintick)
    line([bintick(tt) bintick(tt)], [-2 4], 'color','k','LineStyle',':')
end
xlabel('time from go cue(s)')
ylabel('CD projection(a.u.)')

subplot(3,3,8)
hold on
plot(across_CD_response{1,1}','color',[1 0 0 .1],'LineWidth',.1)
hold on
plot(across_CD_response{1,2}','color',[0 0 1 .1],'LineWidth',.1)
hold on
plot(mean(across_CD_response{1,1},1),'color',[1 0 0],'LineWidth',2)
hold on
plot(mean(across_CD_response{1,2},1),'color',[0 0 1],'LineWidth',2)
ylim([-2 4])

xlim(xlimbd)
xticks(bintick)
xticklabels(xlabeltext)
for tt=1:length(bintick)
    line([bintick(tt) bintick(tt)], [-2 4], 'color','k','LineStyle',':')
end
xlabel('time from go cue(s)')
ylabel('CD projection(a.u.)')

subplot(3,3,9)
hold on
plot(across_CD_response{1,3}','color',[1 0 0 .1],'LineWidth',.1)
hold on
plot(across_CD_response{1,4}','color',[0 0 1 .1],'LineWidth',.1)
hold on
plot(mean(across_CD_response{1,3},1),'color',[1 0 0],'LineWidth',2)
hold on
plot(mean(across_CD_response{1,4},1),'color',[0 0 1],'LineWidth',2)
ylim([-2 4])

xlim(xlimbd)
xticks(bintick)
xticklabels(xlabeltext)
for tt=1:length(bintick)
    line([bintick(tt) bintick(tt)], [-2 4], 'color','k','LineStyle',':')
end
xlabel('time from go cue(s)')
ylabel('CD projection(a.u.)')

hold on
sgtitle(strcat('within session decoder accuracy >',num2str(decoder_cri)))

sz = 4;
figure

temp1 = mean(within_decoder_s(:,1));
temp2 = mean(within_decoder_d(:,1));
temp3 = mean(within_decoder_r(:,1));
temp11 = std(within_decoder_s(:,1))/sqrt(length(within_decoder_s));
temp12 = std(within_decoder_d(:,1))/sqrt(length(within_decoder_d));
temp13 = std(within_decoder_r(:,1))/sqrt(length(within_decoder_r));

hold on
scatter(rand(1,s1),within_decoder_s(:,1),sz,'k')
hold on
scatter(rand(1,d1)+7,within_decoder_d(:,1),sz,'k')
hold on
scatter(rand(1,r1)+14,within_decoder_r(:,1),sz,'k')
hold on
line([1 2], [temp1 temp1], 'color', 'k')
line([1.5 1.5], [temp1-temp11 temp1+temp11], 'color', 'k') 
line([8 9], [temp2 temp2], 'color', 'k')
line([8.5 8.5], [temp2-temp12 temp2+temp12], 'color', 'k') 
line([15 16], [temp3 temp3], 'color', 'k')
line([15.5 15.5], [temp3-temp13 temp3+temp13], 'color', 'k') 

temp1 = mean(across_decoder_s(:,1));
temp2 = mean(across_decoder_d(:,1));
temp3 = mean(across_decoder_r(:,1));
temp11 = std(across_decoder_s(:,1))/sqrt(length(across_decoder_s));
temp12 = std(across_decoder_d(:,1))/sqrt(length(across_decoder_d));
temp13 = std(across_decoder_r(:,1))/sqrt(length(across_decoder_r));

hold on
scatter(rand(1,length(across_decoder_s))+2,across_decoder_s(:,1),sz,'b')
hold on
scatter(rand(1,length(across_decoder_d))+9,across_decoder_d(:,1),sz,'b')
hold on
scatter(rand(1,length(across_decoder_r))+16,across_decoder_r(:,1),sz,'b')
hold on
line([3 4], [temp1 temp1], 'color', 'b')
line([3.5 3.5], [temp1-temp11 temp1+temp11], 'color', 'b') 
line([10 11], [temp2 temp2], 'color', 'b')
line([10.5 10.5], [temp2-temp12 temp2+temp12], 'color', 'b') 
line([17 18], [temp3 temp3], 'color', 'b')
line([17.5 17.5], [temp3-temp13 temp3+temp13], 'color', 'b') 

temp1 = mean(across_decoder_s(:,2));
temp2 = mean(across_decoder_d(:,2));
temp3 = mean(across_decoder_r(:,2));
temp11 = std(across_decoder_s(:,2))/sqrt(length(across_decoder_s));
temp12 = std(across_decoder_d(:,2))/sqrt(length(across_decoder_d));
temp13 = std(across_decoder_r(:,2))/sqrt(length(across_decoder_r));

hold on
scatter(rand(1,length(across_decoder_s))+4,across_decoder_s(:,2),sz,'b')
hold on
scatter(rand(1,length(across_decoder_d))+11,across_decoder_d(:,2),sz,'b')
hold on
scatter(rand(1,length(across_decoder_r))+18,across_decoder_r(:,2),sz,'b')
hold on
line([5 6], [temp1 temp1], 'color', 'b')
line([5.5 5.5], [temp1-temp11 temp1+temp11], 'color', 'b') 
line([12 13], [temp2 temp2], 'color', 'b')
line([12.5 12.5], [temp2-temp12 temp2+temp12], 'color', 'b') 
line([19 20], [temp3 temp3], 'color', 'b')
line([19.5 19.5], [temp3-temp13 temp3+temp13], 'color', 'b') 

xticklabeltext={'S-w/i','S-LE','S-EL','D-w/i','D-LE','D-EL','R-w/i','R-LE','R-EL'};
ylim([0 1]);
xlim([-2 22])
ylabel('accuracy')
xticks([1 3 5 8 10 12 15 17 19]) 
xticklabels(xticklabeltext)

hold on
sgtitle(strcat('within session decoder accuracy >',num2str(decoder_cri)))

sample_bin = round([1.57 2.87]*6);
delay_bin = round([2.87 4.17]*6);
response_bin = round([4.27 6.27]*6);
    
for i=1:size(within_CD_sample{1,1},1)
    delta_CD_sample_within(i,1) = mean(within_CD_sample{1,2}(i,sample_bin(1):sample_bin(2)))-mean(within_CD_sample{1,1}(i,sample_bin(1):sample_bin(2)));
end

for i=1:size(within_CD_delay{1,1},1)
    delta_CD_delay_within(i,1) = mean(within_CD_delay{1,2}(i,delay_bin(1):delay_bin(2)))-mean(within_CD_delay{1,1}(i,delay_bin(1):delay_bin(2)));
end

for i=1:size(within_CD_response{1,1},1)
    delta_CD_response_within(i,1) = mean(within_CD_response{1,2}(i,response_bin(1):response_bin(2)))-mean(within_CD_response{1,1}(i,response_bin(1):response_bin(2)));
end

for i=1:size(across_CD_sample{1,1},1)
    delta_CD_sample_across1(i,1) = mean(across_CD_sample{1,2}(i,sample_bin(1):sample_bin(2)))-mean(across_CD_sample{1,1}(i,sample_bin(1):sample_bin(2)));
end

for i=1:size(across_CD_delay{1,1},1)
    delta_CD_delay_across1(i,1) = mean(across_CD_delay{1,2}(i,delay_bin(1):delay_bin(2)))-mean(across_CD_delay{1,1}(i,delay_bin(1):delay_bin(2)));
end

for i=1:size(across_CD_response{1,1},1)
    delta_CD_response_across1(i,1) = mean(across_CD_response{1,2}(i,response_bin(1):response_bin(2)))-mean(across_CD_response{1,1}(i,response_bin(1):response_bin(2)));
end

for i=1:size(across_CD_sample{1,1},1)
    delta_CD_sample_across2(i,1) = mean(across_CD_sample{1,4}(i,sample_bin(1):sample_bin(2)))-mean(across_CD_sample{1,3}(i,sample_bin(1):sample_bin(2)));
end

for i=1:size(across_CD_delay{1,1},1)
    delta_CD_delay_across2(i,1) = mean(across_CD_delay{1,4}(i,delay_bin(1):delay_bin(2)))-mean(across_CD_delay{1,3}(i,delay_bin(1):delay_bin(2)));
end

for i=1:size(across_CD_response{1,1},1)
    delta_CD_response_across2(i,1) = mean(across_CD_response{1,4}(i,response_bin(1):response_bin(2)))-mean(across_CD_response{1,3}(i,response_bin(1):response_bin(2)));
end

figure
subplot(1,3,1)
hold on
scatter(rand(length(delta_CD_sample_within),1),delta_CD_sample_within,sz,'k')
hold on
scatter(rand(length(delta_CD_sample_across1),1)+2,delta_CD_sample_across1,sz,'b')
hold on
scatter(rand(length(delta_CD_sample_across2),1)+4,delta_CD_sample_across2,sz,'b')

clear temp1 temp2 temp3 temp11 temp12 temp13
temp1 = mean(delta_CD_sample_within);
temp2 = mean(delta_CD_sample_across1);
temp3 = mean(delta_CD_sample_across2);
temp11 = std(delta_CD_sample_within)/sqrt(length(delta_CD_sample_within));
temp12 = std(delta_CD_sample_across1)/sqrt(length(delta_CD_sample_across1));
temp13 = std(delta_CD_sample_across2)/sqrt(length(delta_CD_sample_across2));
hold on
line([1 2], [temp1 temp1], 'color', 'k')
line([1.5 1.5], [temp1-temp11 temp1+temp11], 'color', 'k') 
line([3 4], [temp2 temp2], 'color', 'b')
line([3.5 3.5], [temp2-temp12 temp2+temp12], 'color', 'b') 
line([5 6], [temp3 temp3], 'color', 'b')
line([5.5 5.5], [temp3-temp13 temp3+temp13], 'color', 'b') 
xlim([-1 7])
ylabel('delta CD projection')
xticks([1 3 5]) 
xticklabels({'S-w/i','S-LE','S-EL'})
hold on
title('sample CD')

subplot(1,3,2)
hold on
scatter(rand(length(delta_CD_delay_within),1),delta_CD_delay_within,sz,'k')
hold on
scatter(rand(length(delta_CD_delay_across1),1)+2,delta_CD_delay_across1,sz,'b')
hold on
scatter(rand(length(delta_CD_delay_across2),1)+4,delta_CD_delay_across2,sz,'b')

clear temp1 temp2 temp3 temp11 temp12 temp13
temp1 = mean(delta_CD_delay_within);
temp2 = mean(delta_CD_delay_across1);
temp3 = mean(delta_CD_delay_across2);
temp11 = std(delta_CD_delay_within)/sqrt(length(delta_CD_delay_within));
temp12 = std(delta_CD_delay_across1)/sqrt(length(delta_CD_delay_across1));
temp13 = std(delta_CD_delay_across2)/sqrt(length(delta_CD_delay_across2));
hold on
line([1 2], [temp1 temp1], 'color', 'k')
line([1.5 1.5], [temp1-temp11 temp1+temp11], 'color', 'k') 
line([3 4], [temp2 temp2], 'color', 'b')
line([3.5 3.5], [temp2-temp12 temp2+temp12], 'color', 'b') 
line([5 6], [temp3 temp3], 'color', 'b')
line([5.5 5.5], [temp3-temp13 temp3+temp13], 'color', 'b') 

xlim([-1 7])
ylabel('delta CD projection')
xticks([1 3 5]) 
xticklabels({'D-w/i','D-LE','D-EL'})
hold on
title('delay CD')

subplot(1,3,3)
hold on
scatter(rand(length(delta_CD_response_within),1),delta_CD_response_within,sz,'k')
hold on
scatter(rand(length(delta_CD_response_across1),1)+2,delta_CD_response_across1,sz,'b')
hold on
scatter(rand(length(delta_CD_response_across2),1)+4,delta_CD_response_across2,sz,'b')

clear temp1 temp2 temp3 temp11 temp12 temp13
temp1 = mean(delta_CD_response_within);
temp2 = mean(delta_CD_response_across1);
temp3 = mean(delta_CD_response_across2);
temp11 = std(delta_CD_response_within)/sqrt(length(delta_CD_response_within));
temp12 = std(delta_CD_response_across1)/sqrt(length(delta_CD_response_across1));
temp13 = std(delta_CD_response_across2)/sqrt(length(delta_CD_response_across2));
hold on
line([1 2], [temp1 temp1], 'color', 'k')
line([1.5 1.5], [temp1-temp11 temp1+temp11], 'color', 'k') 
line([3 4], [temp2 temp2], 'color', 'b')
line([3.5 3.5], [temp2-temp12 temp2+temp12], 'color', 'b') 
line([5 6], [temp3 temp3], 'color', 'b')
line([5.5 5.5], [temp3-temp13 temp3+temp13], 'color', 'b') 

xlim([-1 7])
ylabel('delta CD projection')
xticks([1 3 5]) 
xticklabels({'R-w/i','R-LE','R-EL'})
hold on
title('response CD')


sz = 5;
figure
subplot(2,2,1)
hold on
scatter(rand(1,length(within_psth_sim))*-5,within_psth_sim(:,1),sz,'k')
hold on
scatter(across_psth_sim(:,3),across_psth_sim(:,1),sz,'b')
ylim([0 1])
xlim([-8 48])
xlabel('delta days')
ylabel('corr coef (r)')
title('PSTH similiarty - AND gate')

subplot(2,2,2)
hold on
scatter(rand(1,length(within_psth_sim))*-5,within_psth_sim(:,2),sz,'k')
hold on
scatter(across_psth_sim(:,3),across_psth_sim(:,2),sz,'b')
ylim([0 1])
xlim([-8 48])
xlabel('delta days')
ylabel('corr coef (r)')
title('PSTH similiarty - NAND gate')

subplot(2,2,3)
hold on
scatter(rand(1,length(within_popsel_sim))*-5,within_popsel_sim(:,1),sz,'k')
hold on
scatter(across_popsel_sim(:,3),across_popsel_sim(:,1),sz,'b')
ylim([0 1])
xlim([-8 48])
xlabel('delta days')
ylabel('corr coef (r)')
title('pop sel vec similiarty - AND gate')

subplot(2,2,4)
hold on
scatter(rand(1,length(within_popsel_sim))*-5,within_popsel_sim(:,2),sz,'k')
hold on
scatter(across_popsel_sim(:,3),across_popsel_sim(:,2),sz,'b')
ylim([0 1])
xlim([-8 48])
xlabel('delta days')
ylabel('corr coef (r)')
title('pop sel vec similiarty - NAND gate')