
binsize = 0.05;
edges1 = [-1-binsize/2:binsize:1+binsize/2];

figure

subplot(2,3,1)
clear targetdata metrics_matrix2 metrics_matrix
targetdata = pop_sel_s_re{1,2};

for i=1:size(edges1,2)-1         
    for j=1:size(edges1,2)-1
        metrics_matrix(i,j) = length(find(targetdata(:,1) > edges1(i) & targetdata(:,1) <= edges1(i+1) & targetdata(:,2) > edges1(j) & targetdata(:,2) <= edges1(j+1)));
    end
end
metrics_matrix2 = flip(metrics_matrix,2)';
imagesc(metrics_matrix2)
xlabel('Weight(early)','fontsize',12)
ylabel('Weight(late)','fontsize',12)
hold on
title('sample')

subplot(2,3,2)
clear targetdata
targetdata = pop_sel_d_re{1,3};
for i=1:size(edges1,2)-1         
    for j=1:size(edges1,2)-1
        metrics_matrix(i,j) = length(find(targetdata(:,1) > edges1(i) & targetdata(:,1) <= edges1(i+1) & targetdata(:,2) > edges1(j) & targetdata(:,2) <= edges1(j+1)));
    end
end
metrics_matrix2 = flip(metrics_matrix,2)';
imagesc(metrics_matrix2)
xlabel('Weight(early)','fontsize',12)
ylabel('Weight(late)','fontsize',12)
hold on
title('delay')

subplot(2,3,3)
clear targetdata
targetdata = pop_sel_r_re{1,4};
for i=1:size(edges1,2)-1         
    for j=1:size(edges1,2)-1
        metrics_matrix(i,j) = length(find(targetdata(:,1) > edges1(i) & targetdata(:,1) <= edges1(i+1) & targetdata(:,2) > edges1(j) & targetdata(:,2) <= edges1(j+1)));
    end
end
metrics_matrix2 = flip(metrics_matrix,2)';
imagesc(metrics_matrix2)
xlabel('Weight(early)','fontsize',12)
ylabel('Weight(late)','fontsize',12)
hold on
title('response')

subplot(2,3,4)
clear targetdata metrics_matrix2 metrics_matrix
targetdata = resp_all_s;
for i=1:size(edges1,2)-1         
    for j=1:size(edges1,2)-1
        metrics_matrix(i,j) = length(find(targetdata(:,1) > edges1(i) & targetdata(:,1) <= edges1(i+1) & targetdata(:,2) > edges1(j) & targetdata(:,2) <= edges1(j+1)));
    end
end
metrics_matrix2 = flip(metrics_matrix,2)';
imagesc(metrics_matrix2)
xlabel('Weight(early)','fontsize',12)
ylabel('Weight(late)','fontsize',12)
hold on
title('sample')

subplot(2,3,5)
clear targetdata
targetdata = resp_all_d;
for i=1:size(edges1,2)-1         
    for j=1:size(edges1,2)-1
        metrics_matrix(i,j) = length(find(targetdata(:,1) > edges1(i) & targetdata(:,1) <= edges1(i+1) & targetdata(:,2) > edges1(j) & targetdata(:,2) <= edges1(j+1)));
    end
end
metrics_matrix2 = flip(metrics_matrix,2)';
imagesc(metrics_matrix2)
xlabel('Weight(early)','fontsize',12)
ylabel('Weight(late)','fontsize',12)
hold on
title('delay')

subplot(2,3,6)
clear targetdata
targetdata = resp_all_r;
for i=1:size(edges1,2)-1         
    for j=1:size(edges1,2)-1
        metrics_matrix(i,j) = length(find(targetdata(:,1) > edges1(i) & targetdata(:,1) <= edges1(i+1) & targetdata(:,2) > edges1(j) & targetdata(:,2) <= edges1(j+1)));
    end
end
metrics_matrix2 = flip(metrics_matrix,2)';
imagesc(metrics_matrix2)
xlabel('Weight(early)','fontsize',12)
ylabel('Weight(late)','fontsize',12)
hold on
title('response')

colormap(jet)