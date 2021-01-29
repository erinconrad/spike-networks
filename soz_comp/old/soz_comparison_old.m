function metrics = soz_comparison(metrics,is_spike_soz,met,out_folder)

%% Get an aggregate of all spikes based on whether they are in soz or not
nfreq = length(metrics.time.freq);
for f = 1:nfreq
    curr_met = metrics.time.freq(f).(met);
    
    in = [];
    out = [];
    
    for p = 1:length(curr_met.pt)
        zs = curr_met.pt(p).spike.zs;
        is_soz = logical(is_spike_soz(p).is_soz);
        soz_zs = zs(is_soz,:);
        not_zs = zs(~is_soz,:);
        in = [in;soz_zs];
        out = [out;not_zs];
    end
    
    metrics.time.freq(f).(met).in = in;
    metrics.time.freq(f).(met).out = out;
    
    % Compare the last time point for which the majority are not nans
    num_nans = sum(isnan(in),1) + sum(isnan(out),1);
    perc_nans = num_nans./repmat(size(in,1)+size(out,1),1,size(in,2));
    minority_nans = perc_nans < 0.5;
    last_minority_nan = find(minority_nans);
    last_minority_nan = last_minority_nan(end);
    
    % Do an independent two sample t test comparing the zs at that point
    [~,pval,~,stats] = ttest2(in(:,last_minority_nan),out(:,last_minority_nan));
    
    metrics.time.freq(f).(met).pval = pval;
    metrics.time.freq(f).(met).stats = stats;
    metrics.time.freq(f).(met).minority_nans = minority_nans;
end

%% Plot
all_p = cell((nfreq)*2,2);
pcount = 0;
z_range = [0 0];
figure
set(gcf,'position',[1 100 1450 300])
[ha, pos] = tight_subplot(1, nfreq, [0.15 0.01], [0.2 0.11], [0.05 0.01]);
for f = 1:nfreq
    curr_met = metrics.time.freq(f).(met);
    in = curr_met.in;
    out = curr_met.out;
    times = metrics.time.freq(1).sd.pt(1).times;
    minority_nans = metrics.time.freq(f).(met).minority_nans;
    
    % Only take minority nan times
    in = in(:,minority_nans);
    out = out(:,minority_nans);
    times = times(minority_nans);
    
    in_mean = nanmean(in,1);
    out_mean = nanmean(out,1);
    in_se = nanstd(in,0,1)./sqrt(sum(~isnan(in),1));
    out_se = nanstd(out,0,1)./sqrt(sum(~isnan(out),1));
    pval = curr_met.pval;
    
    axes(ha(f))
    errorbar(times,in_mean,in_se,'linewidth',2);
    hold on
    errorbar(times,out_mean,out_se,'linewidth',2);
    
    maxyloc = [in_mean+in_se,...
        out_mean+out_se,...
    in_mean-in_se,...
    out_mean-out_se];
    
    % adjust z_range
    if max(maxyloc) > z_range(2)
        z_range(2) = max(maxyloc);
    end

    if min(maxyloc) < z_range(1)
        z_range(1) = min(maxyloc);
    end
    
  
    pcount = pcount + 1;
    all_p{pcount,2} = f;
    all_p{pcount,1} = pval;
    
    if f == 3
        legend({'Seizure onset','Not seizure onset'},'fontsize',20)
    end
        
end

for f = 1:nfreq
    axes(ha(f))
    ylim(z_range)
    
    if f ~= 1
        yticklabels([])
    elseif f == 1
        ylabel('Normalized power')
    end
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    prettyp = pretty_p(all_p{f,1},3);
    text((xl(1)+xl(2))/2,yl(1)+0.9*(yl(2)-yl(1)),...
        sprintf('%s',prettyp),'fontsize',20,...
        'horizontalalignment','center')
    
    if f == 1
        title('Sub-gamma (<30 Hz)');
    elseif f == 2
        title('Low gamma (30-100 Hz)');
        xlabel('Time relative to spike peak (s)');
    elseif f == 3
        title('High gamma (>100 Hz)');
    end
    
    set(gca,'fontsize',20)
end

print(gcf,[out_folder,'soz_comparison'],'-depsc');

end