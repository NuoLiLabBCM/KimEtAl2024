%%
sz = 20;
colorcode = hsv(size(decoder_s_re,1));
figure
for i=1:size(colorcode,1)
    subplot(3,1,1)
    hold on
    scatter(rand(size(decoder_s_re{i,6},1),1)*-5,decoder_s_re{i,6},sz,colorcode(i,:),'filled','MarkerFaceAlpha',.2)
    hold on
    scatter(decoder_s_re{i,1},decoder_s_re{i,7},sz,colorcode(i,:),'filled','MarkerFaceAlpha',.7)

    subplot(3,1,2)
    hold on
    scatter(rand(size(decoder_d_re{i,6},1),1)*-5,decoder_d_re{i,6},sz,colorcode(i,:),'filled','MarkerFaceAlpha',.2)
    hold on
    scatter(decoder_d_re{i,1},decoder_d_re{i,7},sz,colorcode(i,:),'filled','MarkerFaceAlpha',.7)
    
    subplot(3,1,3)
    hold on
    scatter(rand(size(decoder_r_re{i,6},1),1)*-5,decoder_r_re{i,6},sz,colorcode(i,:),'filled','MarkerFaceAlpha',.2)
    hold on
    scatter(decoder_r_re{i,1},decoder_r_re{i,7},sz,colorcode(i,:),'filled','MarkerFaceAlpha',.7)
    
    
end

for i=1:3
    subplot(3,1,i)
    ylim([0 1])
    ylabel('Decoding accuracy','fontsize',15)
    xlim([-8 62])
    xlabel('Delta days','fontsize',15)
    hold on
    line([0 62], [.5 .5],'color','k','LineStyle',':')
    ax=gca;
    ax.XAxis.FontSize=11;
    ax.YAxis.FontSize=11;
end


subplot(3,1,1)
hold on
title(strcat('CD Sample, 8 mice,',num2str(round(length(within_decoder_s)/2)),' FOVs'),'fontsize',12)

subplot(3,1,2)
hold on
title(strcat('CD Delay, 8 mice,',num2str(round(length(within_decoder_d)/2)),' FOVs'),'fontsize',12)

subplot(3,1,3)
hold on
title(strcat('CD Response, 8 mice,',num2str(round(length(within_decoder_r)/2)),' FOVs'),'fontsize',12)

set(gcf,'color','w')

hold on
sgtitle('averaging two CD decoder')
%exportgraphics(gcf,'SameContextDAovertime_SDR.emf','ContentType','vector')