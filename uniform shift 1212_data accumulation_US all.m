clear; clc;

ep = 0;

default_cd = cd;

cd('uniform shift 1212_US all')
disp('data accumulation')
currdir = dir;
for i=1:size(currdir,1)
    
    currfn = currdir(i).name;
    
    if length(find(currfn == '_')) == 5
        
        load(currfn)
        disp(currfn)
        ep = ep + 1;
        end_points_all{ep,1} = end_points;
        vectors_all{ep,1} = US_1;
        vectors_all{ep,2} = US_2;
        vectors_all{ep,3} = US_3;
        vectors_all{ep,4} = CD_1;
        vectors_all{ep,5} = CD_2;
        vectors_all{ep,6} = CD_3;
        vectors_all{ep,7} = CD_4;
        vectors_dot(ep,1) = corr(CD_1,CD_3);
        vectors_dot(ep,2) = corr(CD_2,CD_4);
        vectors_dot(ep,3) = corr(US_1,US_2);
        vectors_dot(ep,4) = corr(US_2,US_3);
        vectors_dot(ep,5) = corr(US_1,US_3);
        %vectors_dot(ep,5) = corr(CD_1,CD_2);        
        
        fn_al{ep,1} = fn2;
        mouseid(ep,1) = str2num(fn2(6:8));
        clear end_points
    end
end

cd(default_cd)
disp('data saving')
save T1T2T1T2_US_all.mat end_points_all fn_al vectors_dot vectors_all mouseid

disp('dot product plotting')
mouseidpool = [122 130 134];
colorcode = hsv(length(mouseidpool));

vecter_label = {'CD1 vs. CD1r','CD2 vs. CD2r','US1 vs. US2','US1 vs. US3'};

figure
for i=1:length(mouseidpool)    
    hold on
    plot(vectors_dot(find(mouseid == mouseidpool(i)),[1 2 3 5])','color',colorcode(i,:),'LineWidth',1)
    
end
ylim([-1 1])
xlim([0.5 4.5])
xticks([1:1:4])
xticklabels(vecter_label)
line([-.5 4.5],[0 0],'color','k','LineStyle',':')
tempcri = [1 2 3 5];
for i=1:4
    clear temp1 temp2
    temp1 = mean(vectors_dot(:,tempcri(i)));
    temp2 = std(vectors_dot(:,tempcri(i)))/sqrt(length(mouseid));
    
    hold on
    line([i i],[temp1+temp2 temp1-temp2],'color','k','LineWidth',2)
    line([i-.2 i+.2],[temp1 temp1],'color','k','LineWidth',2)
    
end
hold on
sgtitle('dot products between vectors')

vecter_label = {'CD1 vs. CD1r','CD2 vs. CD2r','US1 vs. US2','US2 vs. US3','US1 vs. US3'};

figure
for i=1:length(mouseidpool)    
    hold on
    plot(vectors_dot(find(mouseid == mouseidpool(i)),:)','color',colorcode(i,:),'LineWidth',1)
    
end
ylim([-1 1])
xlim([0.5 5.5])
xticks([1:1:5])
xticklabels(vecter_label)
line([-.5 5.5],[0 0],'color','k','LineStyle',':')

for i=1:5
    clear temp1 temp2
    temp1 = mean(vectors_dot(:,i));
    temp2 = std(vectors_dot(:,i))/sqrt(length(mouseid));
    
    hold on
    line([i i],[temp1+temp2 temp1-temp2],'color','k','LineWidth',2)
    line([i-.2 i+.2],[temp1 temp1],'color','k','LineWidth',2)
    
end
hold on
sgtitle('dot products between vectors')
%%




%%
disp('end point data processing')
for i=1:size(end_points_all,1)
    
    %%% Trial averaged original end points
    % T1 dataset
    % lick right
    data_combined{1,1}(i,1) = mean(end_points_all{i,1}{1,1}{7,1}(:,1));
    data_combined{1,1}(i,2) = mean(end_points_all{i,1}{1,1}{7,1}(:,2));
    % lick left
    data_combined{1,2}(i,1) = mean(end_points_all{i,1}{1,2}{7,1}(:,1));
    data_combined{1,2}(i,2) = mean(end_points_all{i,1}{1,2}{7,1}(:,2));
    
    % T2 dataset
    % lick right
    data_combined{1,3}(i,1) = mean(end_points_all{i,1}{2,1}{7,1}(:,1));
    data_combined{1,3}(i,2) = mean(end_points_all{i,1}{2,1}{7,1}(:,2));
    % lick left
    data_combined{1,4}(i,1) = mean(end_points_all{i,1}{2,2}{7,1}(:,1));
    data_combined{1,4}(i,2) = mean(end_points_all{i,1}{2,2}{7,1}(:,2));
    
    % T1' dataset
    % lick right
    data_combined{1,5}(i,1) = mean(end_points_all{i,1}{3,1}{7,1}(:,1));
    data_combined{1,5}(i,2) = mean(end_points_all{i,1}{3,1}{7,1}(:,2));
    % lick left
    data_combined{1,6}(i,1) = mean(end_points_all{i,1}{3,2}{7,1}(:,1));
    data_combined{1,6}(i,2) = mean(end_points_all{i,1}{3,2}{7,1}(:,2));
    
    % T2' dataset
    % lick right
    data_combined{1,7}(i,1) = mean(end_points_all{i,1}{4,1}{7,1}(:,1));
    data_combined{1,7}(i,2) = mean(end_points_all{i,1}{4,1}{7,1}(:,2));
    % lick left
    data_combined{1,8}(i,1) = mean(end_points_all{i,1}{4,2}{7,1}(:,1));
    data_combined{1,8}(i,2) = mean(end_points_all{i,1}{4,2}{7,1}(:,2));
end

% x axis
right_combined{1,1} = horzcat(data_combined{1,1}(:,1),data_combined{1,3}(:,1),data_combined{1,5}(:,1),data_combined{1,7}(:,1));
% y axis
right_combined{1,2} = horzcat(data_combined{1,1}(:,2),data_combined{1,3}(:,2),data_combined{1,5}(:,2),data_combined{1,7}(:,2));

% x axis
left_combined{1,1} = horzcat(data_combined{1,2}(:,1),data_combined{1,4}(:,1),data_combined{1,6}(:,1),data_combined{1,8}(:,1));
% y axis
left_combined{1,2} = horzcat(data_combined{1,2}(:,2),data_combined{1,4}(:,2),data_combined{1,6}(:,2),data_combined{1,8}(:,2));


sz = 50;
figure
subplot(2,5,1)
hold on
% lick right (T1)
scatter(data_combined{1,1}(:,1),data_combined{1,1}(:,2),sz,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerFaceAlpha',.2)

% lick left (T1)
scatter(data_combined{1,2}(:,1),data_combined{1,2}(:,2),sz,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerFaceAlpha',.2)

line([-3 3],[0 0],'color','k','LineStyle',':')
line([0 0],[-3 3],'color','k','LineStyle',':')
xlim([-1.5 2])
ylim([-.9 1.6])
xlabel('CD1')
ylabel('CD2')
hold on
title('all FOVs, T1')

subplot(2,5,2)
hold on
% lick right (T2)
scatter(data_combined{1,3}(:,1),data_combined{1,3}(:,2),sz,'MarkerEdgeColor','c','MarkerFaceColor','c','MarkerFaceAlpha',.2)

% lick left (T2)
scatter(data_combined{1,4}(:,1),data_combined{1,4}(:,2),sz,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerFaceAlpha',.2)

line([-3 3],[0 0],'color','k','LineStyle',':')
line([0 0],[-3 3],'color','k','LineStyle',':')
xlim([-1.5 2])
ylim([-.9 1.6])
ylabel('CD1(ortho)')
xlabel('US')
hold on
title('all FOVs, T2')

subplot(2,5,6)
hold on
% lick right (T1')
scatter(data_combined{1,5}(:,1),data_combined{1,5}(:,2),sz,'x','MarkerEdgeColor','b','MarkerFaceColor','none')

% lick left (T1')
scatter(data_combined{1,6}(:,1),data_combined{1,6}(:,2),sz,'x','MarkerEdgeColor','r','MarkerFaceColor','none')

line([-3 3],[0 0],'color','k','LineStyle',':')
line([0 0],[-3 3],'color','k','LineStyle',':')
xlim([-1.5 2])
ylim([-.9 1.6])
ylabel('CD1(ortho)')
xlabel('US')
hold on
title('all FOVs, T1r')

subplot(2,5,7)
hold on
% lick right (T2)
scatter(data_combined{1,7}(:,1),data_combined{1,7}(:,2),sz,'x','MarkerEdgeColor','c','MarkerFaceColor','none')

% lick left (T2)
scatter(data_combined{1,8}(:,1),data_combined{1,8}(:,2),sz,'x','MarkerEdgeColor','m','MarkerFaceColor','none')

line([-3 3],[0 0],'color','k','LineStyle',':')
line([0 0],[-3 3],'color','k','LineStyle',':')
xlim([-1.5 2])
ylim([-.9 1.6])
ylabel('CD1(ortho)')
xlabel('US')
hold on
title('all FOVs, T2r')

subplot(2,5,3)
hold on
% lick right (T1)
scatter(data_combined{1,1}(:,1),data_combined{1,1}(:,2),sz,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerFaceAlpha',.3)

% lick right (T2)
scatter(data_combined{1,3}(:,1),data_combined{1,3}(:,2),sz,'MarkerEdgeColor','c','MarkerFaceColor','c','MarkerFaceAlpha',.3)

% lick right (T1')
scatter(data_combined{1,5}(:,1),data_combined{1,5}(:,2),sz,'x','MarkerEdgeColor','b','MarkerFaceColor','none')

% lick right (T2')
scatter(data_combined{1,7}(:,1),data_combined{1,7}(:,2),sz,'x','MarkerEdgeColor','c','MarkerFaceColor','none')

line([-3 3],[0 0],'color','k','LineStyle',':')
line([0 0],[-3 3],'color','k','LineStyle',':')
xlim([-1.5 2])
ylim([-.9 1.6])
ylabel('CD1(ortho)')
xlabel('US')
hold on
title('right trials only')

subplot(2,5,8)
hold on
% lick left (T1)
scatter(data_combined{1,2}(:,1),data_combined{1,2}(:,2),sz,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerFaceAlpha',.2)

% lick left (T2)
scatter(data_combined{1,4}(:,1),data_combined{1,4}(:,2),sz,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerFaceAlpha',.2)

% lick left (T1')
scatter(data_combined{1,6}(:,1),data_combined{1,6}(:,2),sz,'x','MarkerEdgeColor','r','MarkerFaceColor','none')

% lick left (T2')
scatter(data_combined{1,8}(:,1),data_combined{1,8}(:,2),sz,'x','MarkerEdgeColor','m','MarkerFaceColor','none')

line([-3 3],[0 0],'color','k','LineStyle',':')
line([0 0],[-3 3],'color','k','LineStyle',':')
xlim([-1.5 2])
ylim([-.9 1.6])
ylabel('CD1(ortho)')
xlabel('US')
hold on
title('left trials only')

xtext_label = {'T1','T2','T1r','T2r'};

subplot(2,5,4)
hold on
plot(right_combined{1,1}','color',[0 0 1 .5],'LineWidth',0.5)
xlim([0.5 4.5])
xticks([1:1:4])
xticklabels(xtext_label)
line([-.5 4.5],[0 0],'color','k','LineStyle',':')

for i=1:4
    clear temp1 temp2
    temp1 = mean(right_combined{1,1}(:,i));
    temp2 = std(right_combined{1,1}(:,i))/sqrt(length(mouseid));
    
    hold on
    line([i i],[temp1+temp2 temp1-temp2],'color','b','LineWidth',1.5)
    line([i-.2 i+.2],[temp1 temp1],'color','b','LineWidth',1.5)
    
end
ylabel('US')
hold on
title('US(x)')

subplot(2,5,5)
hold on
plot(right_combined{1,2}','color',[0 0 1 .5],'LineWidth',0.5)
xlim([0.5 4.5])
xticks([1:1:4])
xticklabels(xtext_label)
line([-.5 4.5],[0 0],'color','k','LineStyle',':')

for i=1:4
    clear temp1 temp2
    temp1 = mean(right_combined{1,2}(:,i));
    temp2 = std(right_combined{1,2}(:,i))/sqrt(length(mouseid));
    
    hold on
    line([i i],[temp1+temp2 temp1-temp2],'color','b','LineWidth',1.5)
    line([i-.2 i+.2],[temp1 temp1],'color','b','LineWidth',1.5)
    
end
xlabel('CD1')
hold on
title('CD1(ortho,y)')


subplot(2,5,9)
hold on
plot(left_combined{1,1}','color',[1 0 0 .5],'LineWidth',0.5)
xlim([0.5 4.5])
xticks([1:1:4])
xticklabels(xtext_label)
line([-.5 4.5],[0 0],'color','k','LineStyle',':')

for i=1:4
    clear temp1 temp2
    temp1 = mean(left_combined{1,1}(:,i));
    temp2 = std(left_combined{1,1}(:,i))/sqrt(length(mouseid));
    
    hold on
    line([i i],[temp1+temp2 temp1-temp2],'color','r','LineWidth',1.5)
    line([i-.2 i+.2],[temp1 temp1],'color','r','LineWidth',1.5)
    
end
ylabel('US')
hold on
title('US(x)')

subplot(2,5,10)
hold on
plot(left_combined{1,2}','color',[1 0 0 .5],'LineWidth',0.5)
xlim([0.5 4.5])
xticks([1:1:4])
xticklabels(xtext_label)
line([-.5 4.5],[0 0],'color','k','LineStyle',':')

for i=1:4
    clear temp1 temp2
    temp1 = mean(left_combined{1,2}(:,i));
    temp2 = std(left_combined{1,2}(:,i))/sqrt(length(mouseid));
    
    hold on
    line([i i],[temp1+temp2 temp1-temp2],'color','r','LineWidth',1.5)
    line([i-.2 i+.2],[temp1 temp1],'color','r','LineWidth',1.5)
    
end
xlabel('CD1')
hold on
title('CD1(ortho,y)')

hold on
sgtitle('CD1 and US: task context 1212')