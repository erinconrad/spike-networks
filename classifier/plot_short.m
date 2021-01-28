function plot_short(metrics,met,f,earliest_rise)

show_all = 0;

%% Initialize figure
figure
set(gcf,'position',[1 100 700 400])


%% Get data
data = metrics.time.freq(f).(met).auc.short.data;
pvals = metrics.time.freq(f).(met).auc.short.ps;
np = size(data,2);
ntimes = size(data,1);
times = metrics.time.freq(f).(met).auc.times;

%% Denote the earliest visual rise
min_rise = min(earliest_rise,[],2); % earliest reviewer notation for each spike
mean_rise_spikes = mean(min_rise); % mean across spikes
before_rise = zeros(ntimes,1);
last_before_rise = 1;
% Round to lowest 0.1 s
for tt = 1:ntimes
    if mean_rise_spikes > times(tt)
        before_rise(tt) = 1;
        last_before_rise = tt;
    end
end


auc_diff = squeeze(data(:,:,1)-data(:,:,2));
mean_auc_diff = mean(auc_diff,2);
std_diff = std(auc_diff,0,2);

if show_all == 1
    errorbar(times(1:end),...
        mean_auc_diff(1:end)...
        ,std_diff(1:end),'ko','markersize',10)
else
    errorbar(times(1:last_before_rise-1),...
        mean_auc_diff(1:last_before_rise-1)...
        ,std_diff(1:last_before_rise-1),'ko','markersize',10)
end
hold on

%% Formatting

xlabel('Time (s)');


ylabel(sprintf('Relative %s %s summed\nacross pre-spike epoch\nSpike-non spike difference',...
    metrics.time.freq(f).name,met));


set(gca,'fontsize',20)



%% Add p-value
sub_alpha = zeros(ntimes,1);
for tt = 1:last_before_rise-1
    num_left = ntimes-tt+1;
    num_sub_alpha = sum(pvals(tt:end) < 0.05/3);
    if num_sub_alpha == num_left
        sub_alpha(tt) = 1;
    end
end



yl = get(gca,'ylim');
yl2 = yl(1) + 1.1*(yl(2)-yl(1));
ylim([yl(1) yl2]);
yl = get(gca,'ylim');
p_yloc = yl(1) + 0.9*(yl(2)-yl(1));
for tt = 1:ntimes
    if sub_alpha(tt) == 1
        text(times(tt),p_yloc,'*','horizontalalignment','center','fontsize',30);
    end
end

endh = plot([times(last_before_rise) times(last_before_rise)],get(gca,'ylim'),'k--','linewidth',2);

legend(endh,'Visual rise','fontsize',20,'location','northwest')



    
    
end

