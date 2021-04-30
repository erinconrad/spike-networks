function compare_pre_to_post(all,m,f,out_folder)

alpha = 0.05;
np = length(all(1).data.freq(1).pt);
all_pt_dat = nan(np,2);
all_pval = nan(np,1);
times = all(m).data.times;

pre_time_max = -0.3;
post_time_min = 1;

for p = 1:np
    %% get data
    % nsp x ntimes x nchs
    data = all(m).data.freq(f).pt(p).sp_or_not(1).data;
    
    %% divide times into pre and post
    pre = data(:,times <= pre_time_max,:);
    post = data(:,times >= post_time_min,:);
    
    %% within patient test
    % average across times and channels
    all_pre = mean(pre,[2 3]);
    all_post = mean(post,[2 3]);
    
    [~,pval] = ttest(all_pre,all_post);
    all_pval(p) = pval;
    
    %% average for aggregate pt test
    % average across all spikes, times, and channels
    avg_pre = median(pre(:));
    avg_post = median(post(:));
    
    all_pt_dat(p,:) = [avg_pre,avg_post];
    
end

% t-test
[~,pval] = ttest(all_pt_dat(:,1),all_pt_dat(:,2));
text_out = get_asterisks(pval,1);

% Plot
if 1
figure
plot([1 2],[all_pt_dat(:,2),all_pt_dat(:,1)],...
    'k--o','linewidth',2,'markersize',10)
hold on
%{
title(sprintf('Spike: %1.2f, Not: %1.2f, p = %1.3f',...
    mean(all_pt_dat(:,1)),mean(all_pt_dat(:,2)),pval));
%}
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Pre-IED','Post-IED'})
yl = ylim;
ylim([yl(1) yl(1) + 1.1*(yl(2)-yl(1))])
yl = ylim;
plot([1 2],[yl(1)+0.9*(yl(2)-yl(1)) yl(1)+0.9*(yl(2)-yl(1))],'k','linewidth',2)
text(1.5,yl(1) + 0.95*(yl(2)-yl(1)),text_out,'fontsize',20,'horizontalalignment','center')


if m == 1
    ylabel('Average power (\muV^2)')
elseif m == 4
    ylabel('Average node strength')
end
set(gca,'fontsize',20)

fprintf('%d of %d had significant difference\n',...
    sum(all_pval<alpha),np);



end