function plot_avg_adj(adj_avg)

adj = adj_avg(2).adj;
ntimes = size(adj,1);
figure
set(gcf,'position',[1 300 1400 200])
ha = tight_subplot(1,6);
count = 0;
for i = 9:14
    count = count + 1;
    axes(ha(count))
    imagesc(squeeze(adj(i,:,:))-squeeze(adj(9,:,:)))
end


end