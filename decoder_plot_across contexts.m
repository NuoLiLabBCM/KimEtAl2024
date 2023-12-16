


%%
sz = 5;
tpcri = 0.15;

figure
subplot(2,3,1)
hold on
scatter(rand(1,length(across_decoder_d(:,1)))/2+0.25,across_decoder_d(:,5),sz,'k')
hold on
scatter(rand(1,length(across_decoder_d(:,2)))/2+1.25,across_decoder_d(:,6),sz,'k')
hold on
scatter(rand(1,length(across_decoder_d(:,1)))/2+2.25,across_decoder_d(:,1),sz,'k')
hold on
scatter(rand(1,length(across_decoder_d(:,2)))/2+3.25,across_decoder_d(:,2),sz,'k')
hold on
clear temp1 temp2
temp1(1,1) = mean(across_decoder_d(:,5));
temp1(1,2) = mean(across_decoder_d(:,6));
temp1(1,3) = mean(across_decoder_d(:,1));
temp1(1,4) = mean(across_decoder_d(:,2));

temp2(1,1) = std(across_decoder_d(:,5))/sqrt(length(across_decoder_d));
temp2(1,2) = std(across_decoder_d(:,6))/sqrt(length(across_decoder_d));
temp2(1,3) = std(across_decoder_d(:,1))/sqrt(length(across_decoder_d));
temp2(1,4) = std(across_decoder_d(:,2))/sqrt(length(across_decoder_d));

for i=1:4
    hold on
    b1 = bar(i-0.5,temp1(1,i),'FaceColor','k');
    b1.FaceAlpha = 0.1;
    hold on
    errorbar(i-0.5,temp1(1,i),temp2(1,i),'color','k','LineWidth',1)
end

ylim([0 1])
ylabel('Decoding accuracy')
xticks([.5:1:3.5])
xticklabels({'11','22','12','21'})
xlim([-.5 4.5])
hold on
title(strcat('Delay CD(n=',num2str(length(across_decoder_d)),')'))

ax=gca;
ax.XAxis.FontSize=11;
ax.YAxis.FontSize=11;

subplot(2,3,4)
hold on
scatter(ones(1,length(across_decoder_d(:,1)))*.7,across_decoder_d(:,5),sz,'k')
hold on
scatter(ones(1,length(across_decoder_d(:,2)))*1.7,across_decoder_d(:,6),sz,'k')
hold on
scatter(ones(1,length(across_decoder_d(:,1)))*2.7,across_decoder_d(:,1),sz,'k')
hold on
scatter(ones(1,length(across_decoder_d(:,2)))*3.7,across_decoder_d(:,2),sz,'k')
for t=1:size(across_decoder_d,1)
    hold on
    line([1 2]-.3, [across_decoder_d(t,5) across_decoder_d(t,6)], 'color',[0 0 0 tpcri])
    line([2 3]-.3, [across_decoder_d(t,6) across_decoder_d(t,1)], 'color',[0 0 0 tpcri])
    line([3 4]-.3, [across_decoder_d(t,1) across_decoder_d(t,2)], 'color',[0 0 0 tpcri])
end
hold on
clear temp1 temp2
temp1(1,1) = mean(across_decoder_d(:,5));
temp1(1,2) = mean(across_decoder_d(:,6));
temp1(1,3) = mean(across_decoder_d(:,1));
temp1(1,4) = mean(across_decoder_d(:,2));

temp2(1,1) = std(across_decoder_d(:,5))/sqrt(length(across_decoder_d));
temp2(1,2) = std(across_decoder_d(:,6))/sqrt(length(across_decoder_d));
temp2(1,3) = std(across_decoder_d(:,1))/sqrt(length(across_decoder_d));
temp2(1,4) = std(across_decoder_d(:,2))/sqrt(length(across_decoder_d));

for i=1:4
    hold on
    b1 = bar(i-0.5,temp1(1,i),'FaceColor','k');
    b1.FaceAlpha = 0.1;
    hold on
    errorbar(i-0.5,temp1(1,i),temp2(1,i),'color','k','LineWidth',1)
end

ylim([0 1])
ylabel('Decoding accuracy')
xticks([.5:1:3.5])
xticklabels({'11','22','12','21'})
xlim([-.5 4.5])
hold on
title(strcat('Delay CD(n=',num2str(length(across_decoder_d)),')'))

ax=gca;
ax.XAxis.FontSize=11;
ax.YAxis.FontSize=11;

subplot(2,3,2)
hold on
scatter(rand(1,length(across_decoder_s(:,1)))/2+0.25,across_decoder_s(:,5),sz,'k')
hold on
scatter(rand(1,length(across_decoder_s(:,2)))/2+1.25,across_decoder_s(:,6),sz,'k')
hold on
scatter(rand(1,length(across_decoder_s(:,1)))/2+2.25,across_decoder_s(:,1),sz,'k')
hold on
scatter(rand(1,length(across_decoder_s(:,2)))/2+3.25,across_decoder_s(:,2),sz,'k')
clear temp1 temp2
temp1(1,1) = mean(across_decoder_s(:,5));
temp1(1,2) = mean(across_decoder_s(:,6));
temp1(1,3) = mean(across_decoder_s(:,1));
temp1(1,4) = mean(across_decoder_s(:,2));

temp2(1,1) = std(across_decoder_s(:,5))/sqrt(length(across_decoder_s));
temp2(1,2) = std(across_decoder_s(:,6))/sqrt(length(across_decoder_s));
temp2(1,3) = std(across_decoder_s(:,1))/sqrt(length(across_decoder_s));
temp2(1,4) = std(across_decoder_s(:,2))/sqrt(length(across_decoder_s));

for i=1:4
    hold on
    b1 = bar(i-0.5,temp1(1,i),'FaceColor','k');
    b1.FaceAlpha = 0.1;
    hold on
    errorbar(i-0.5,temp1(1,i),temp2(1,i),'color','k','LineWidth',1)
end

ylim([0 1])
ylabel('Decoding accuracy')
xticks([.5:1:3.5])
xticklabels({'11','22','12','21'})
xlim([-.5 4.5])
hold on
title(strcat('Sample CD(n=',num2str(length(across_decoder_s)),')'))

ylim([0 1])
ylabel('Decoding accuracy')

ax=gca;
ax.XAxis.FontSize=11;
ax.YAxis.FontSize=11;

subplot(2,3,5)
hold on
scatter(ones(1,length(across_decoder_s(:,1)))*.7,across_decoder_s(:,5),sz,'k')
hold on
scatter(ones(1,length(across_decoder_s(:,2)))*1.7,across_decoder_s(:,6),sz,'k')
hold on
scatter(ones(1,length(across_decoder_s(:,1)))*2.7,across_decoder_s(:,1),sz,'k')
hold on
scatter(ones(1,length(across_decoder_s(:,2)))*3.7,across_decoder_s(:,2),sz,'k')
for t=1:size(across_decoder_s,1)
    hold on
    line([1 2]-.3, [across_decoder_s(t,5) across_decoder_s(t,6)], 'color',[0 0 0 tpcri])
    line([2 3]-.3, [across_decoder_s(t,6) across_decoder_s(t,1)], 'color',[0 0 0 tpcri])
    line([3 4]-.3, [across_decoder_s(t,1) across_decoder_s(t,2)], 'color',[0 0 0 tpcri])
end
hold on
clear temp1 temp2
temp1(1,1) = mean(across_decoder_s(:,5));
temp1(1,2) = mean(across_decoder_s(:,6));
temp1(1,3) = mean(across_decoder_s(:,1));
temp1(1,4) = mean(across_decoder_s(:,2));

temp2(1,1) = std(across_decoder_s(:,5))/sqrt(length(across_decoder_s));
temp2(1,2) = std(across_decoder_s(:,6))/sqrt(length(across_decoder_s));
temp2(1,3) = std(across_decoder_s(:,1))/sqrt(length(across_decoder_s));
temp2(1,4) = std(across_decoder_s(:,2))/sqrt(length(across_decoder_s));

for i=1:4
    hold on
    b1 = bar(i-0.5,temp1(1,i),'FaceColor','k');
    b1.FaceAlpha = 0.1;
    hold on
    errorbar(i-0.5,temp1(1,i),temp2(1,i),'color','k','LineWidth',1)
end

ylim([0 1])
ylabel('Decoding accuracy')
xticks([.5:1:3.5])
xticklabels({'11','22','12','21'})
xlim([-.5 4.5])
hold on
title(strcat('Sample CD(n=',num2str(length(across_decoder_s)),')'))

ax=gca;
ax.XAxis.FontSize=11;
ax.YAxis.FontSize=11;

subplot(2,3,3)
hold on
scatter(rand(1,length(across_decoder_r(:,1)))/2+0.25,across_decoder_r(:,5),sz,'k')
hold on
scatter(rand(1,length(across_decoder_r(:,2)))/2+1.25,across_decoder_r(:,6),sz,'k')
hold on
scatter(rand(1,length(across_decoder_r(:,1)))/2+2.25,across_decoder_r(:,1),sz,'k')
hold on
scatter(rand(1,length(across_decoder_r(:,2)))/2+3.25,across_decoder_r(:,2),sz,'k')
clear temp1 temp2
temp1(1,1) = mean(across_decoder_r(:,5));
temp1(1,2) = mean(across_decoder_r(:,6));
temp1(1,3) = mean(across_decoder_r(:,1));
temp1(1,4) = mean(across_decoder_r(:,2));

temp2(1,1) = std(across_decoder_r(:,5))/sqrt(length(across_decoder_r));
temp2(1,2) = std(across_decoder_r(:,6))/sqrt(length(across_decoder_r));
temp2(1,3) = std(across_decoder_r(:,1))/sqrt(length(across_decoder_r));
temp2(1,4) = std(across_decoder_r(:,2))/sqrt(length(across_decoder_r));

for i=1:4
    hold on
    b1 = bar(i-0.5,temp1(1,i),'FaceColor','k');
    b1.FaceAlpha = 0.1;
    hold on
    errorbar(i-0.5,temp1(1,i),temp2(1,i),'color','k','LineWidth',1)
end

ylim([0 1])
ylabel('Decoding accuracy')
xticks([.5:1:3.5])
xticklabels({'11','22','12','21'})
xlim([-.5 4.5])
hold on
title(strcat('Response CD(n=',num2str(length(across_decoder_r)),')'))

ylim([0 1])
ylabel('Decoding accuracy')

ax=gca;
ax.XAxis.FontSize=11;
ax.YAxis.FontSize=11;

subplot(2,3,6)
hold on
scatter(ones(1,length(across_decoder_r(:,1)))*.7,across_decoder_r(:,5),sz,'k')
hold on
scatter(ones(1,length(across_decoder_r(:,2)))*1.7,across_decoder_r(:,6),sz,'k')
hold on
scatter(ones(1,length(across_decoder_r(:,1)))*2.7,across_decoder_r(:,1),sz,'k')
hold on
scatter(ones(1,length(across_decoder_r(:,2)))*3.7,across_decoder_r(:,2),sz,'k')
for t=1:size(across_decoder_r,1)
    hold on
    line([1 2]-.3, [across_decoder_r(t,5) across_decoder_r(t,6)], 'color',[0 0 0 tpcri])
    line([2 3]-.3, [across_decoder_r(t,6) across_decoder_r(t,1)], 'color',[0 0 0 tpcri])
    line([3 4]-.3, [across_decoder_r(t,1) across_decoder_r(t,2)], 'color',[0 0 0 tpcri])
end
hold on
clear temp1 temp2
temp1(1,1) = mean(across_decoder_r(:,5));
temp1(1,2) = mean(across_decoder_r(:,6));
temp1(1,3) = mean(across_decoder_r(:,1));
temp1(1,4) = mean(across_decoder_r(:,2));

temp2(1,1) = std(across_decoder_r(:,5))/sqrt(length(across_decoder_r));
temp2(1,2) = std(across_decoder_r(:,6))/sqrt(length(across_decoder_r));
temp2(1,3) = std(across_decoder_r(:,1))/sqrt(length(across_decoder_r));
temp2(1,4) = std(across_decoder_r(:,2))/sqrt(length(across_decoder_r));

for i=1:4
    hold on
    b1 = bar(i-0.5,temp1(1,i),'FaceColor','k');
    b1.FaceAlpha = 0.1;
    hold on
    errorbar(i-0.5,temp1(1,i),temp2(1,i),'color','k','LineWidth',1)
end

ylim([0 1])
ylabel('Decoding accuracy')
xticks([.5:1:3.5])
xticklabels({'11','22','12','21'})
xlim([-.5 4.5])
hold on
title(strcat('Response CD(n=',num2str(length(across_decoder_r)),')'))

ax=gca;
ax.XAxis.FontSize=11;
ax.YAxis.FontSize=11;

set(gcf,'color','w')

figure
subplot(1,3,1)
clear temp
temp = horzcat(across_decoder_d(:,5:6),across_decoder_d(:,1:2));
imagesc(temp)
caxis([0 1])
hold on
title('Delay CD')
xticks([.5:1:3.5])
xticklabels({'11','22','12','21'})
ylabel('FOVs')
colorbar

subplot(1,3,2)
clear temp
temp = horzcat(across_decoder_s(:,5:6),across_decoder_s(:,1:2));
imagesc(temp)
caxis([0 1])
hold on
title('Sample CD')
xticks([.5:1:3.5])
xticklabels({'11','22','12','21'})
ylabel('FOVs')
colorbar

subplot(1,3,3)
clear temp
temp = horzcat(across_decoder_r(:,5:6),across_decoder_r(:,1:2));
imagesc(temp)
caxis([0 1])
hold on
title('Response CD')
xticks([.5:1:3.5])
xticklabels({'11','22','12','21'})
ylabel('FOVs')
colormap(gray)
colorbar

set(gcf,'color','w')

hold on
sgtitle('Across contexts')

%%
figure(1)

exportgraphics(gcf,'across_contexts_DA_final2.emf','ContentType','vector')
