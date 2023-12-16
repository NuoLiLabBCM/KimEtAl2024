% mouse age analysis
% aa: N(mouse number) x 16 matrix
% bb: finish day
% cc: surgery day


sz = 50;
figure
for i=1:size(aa,1)    
    hold on
    scatter(cc(i,1),i,sz,'k')    
    clear temp
    temp = max(find(aa(i,:)>0));
    
    for j=1:temp
        if aa(i,j) > 0
            if rem(j,3) == 1
                hold on
                scatter(aa(i,j),i,sz,'k','filled')    
                if aa(i,j+1) > 0
                    line([aa(i,j) aa(i,j+1)],[i,i],'color','k','LineWidth',2)
                else
                    line([aa(i,j) bb(i,1)],[i,i],'color','k','LineWidth',2)
                    hold on
                    scatter(bb(i,1),i,sz,'g','filled')
                end
            elseif rem(j,3) == 2
                hold on
                scatter(aa(i,j),i,sz,'b','filled')     
                if aa(i,j+1) > 0
                    line([aa(i,j) aa(i,j+1)],[i,i],'color','b','LineWidth',2)
                else
                    line([aa(i,j) bb(i,1)],[i,i],'color','b','LineWidth',2)
                    hold on
                    scatter(bb(i,1),i,sz,'g','filled')
                end
            elseif rem(j,3) == 0
                hold on
                scatter(aa(i,j),i,sz,'r','filled')     
                if aa(i,j+1) > 0
                    line([aa(i,j) aa(i,j+1)],[i,i],'color','r','LineWidth',2)
                else
                    line([aa(i,j) bb(i,1)],[i,i],'color','r','LineWidth',2)
                    hold on
                    scatter(bb(i,1),i,sz,'g','filled')
                end
            end
        end
    end
end

mouseid = {'JH35','JH116','JH117','JH118','JH119','JH122','JH123','JH130','JH133','JH134','JH158','JH159','JH160','JH176'};
xlabel('Mouse age (postnatal day)')
ylabel('Mouse ID')
ylim([0 length(mouseid)+1])
yticks([1:1:length(mouseid)])
yticklabels(mouseid)
hold on
title('Experimental timeline across mouse age (N=14)')