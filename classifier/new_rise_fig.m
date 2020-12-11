function new_rise_fig(metrics,met)

nfreq = length(metrics.time.freq);

%% Initialize figure
figure
set(gcf,'position',[1 100 1450 270])
[ha, pos] = tight_subplot(1, nfreq, [0.10 0.01], [0.12 0.11], [0.11 0.01]);

z_range = [0 0];
all_pp = cell(nfreq,1);

for f = 1:nfreq
    axes(ha(f))
    
    %% Get data
    dat_sp = metrics.time.freq(f).(met).median_data.data(:,:,1);
    dat_not = metrics.time.freq(f).(met).median_data.data(:,:,2);
    pval = metrics.time.freq(f).(met).auc.pval;
    times = metrics.time.freq(f).(met).pt(1).times;
    
    
    %% Do the plot
    dat_sp_mean = nanmean(dat_sp,1);
    dat_not_mean = nanmean(dat_not,1);
    dat_sp_se = nanstd(dat_sp,0,1)./sqrt(sum(~isnan(dat_sp),1));
    dat_not_se = nanstd(dat_not,0,1)./sqrt(sum(~isnan(dat_not),1));
    sp_pt = errorbar(times,dat_sp_mean,dat_sp_se,'linewidth',2);
    hold on
    not_pt = errorbar(times,dat_not_mean,dat_not_se,'linewidth',2);
    
    maxyloc = [dat_sp_mean+dat_sp_se,dat_not_mean+dat_not_se,...
    dat_sp_mean-dat_sp_se,dat_not_mean-dat_not_se];

    % adjust z_range
    if max(maxyloc) > z_range(2)
        z_range(2) = max(maxyloc);
    end

    if min(maxyloc) < z_range(1)
        z_range(1) = min(maxyloc);
    end
    
    %% Add the p-value
    pp = pretty_p(pval,nfreq);
    all_pp{f} = pp;
      
end

% Loop through axes and adjust ylims and add p values
for f = 1:nfreq
    axes(ha(f))
    
    yl(1) = z_range(1) - 0.1*(z_range(2)-z_range(1));
    yl(2) = z_range(2) + 0.1*(z_range(2)-z_range(1));
    
    ylim(yl)
    
    xl = get(gca,'xlim');
    
    if f ~= 1
        yticklabels([])
    end
    
    text((xl(1)+xl(2))/2,yl(1)+0.9*(yl(2)-yl(1)),...
        sprintf('%s',all_pp{f}),'fontsize',20,...
        'horizontalalignment','center')
end

end