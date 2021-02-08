function plot_short(metrics,met,f,earliest_rise,out_folder,do_plot)

show_all = 0;

%% Pretty names
if strcmp(met,'sd')
    pretty_name = 'absolute power';
    nfreq = 1;
elseif strcmp(met,'ns_avg')
    pretty_name = 'average node strength';
    nfreq = length(metrics.time.freq);
else
    pretty_name = met;
    nfreq = length(metrics.time.freq);
end


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
        ,std_diff(1:end),'ko','markersize',15,...
        'linewidth',2)
else
    errorbar(times(1:last_before_rise-1),...
        mean_auc_diff(1:last_before_rise-1)...
        ,std_diff(1:last_before_rise-1),'ko','markersize',15,...
        'linewidth',2)
end
hold on

%% Formatting

xlabel('Time (s)');

if strcmp(met,'sd')
    ylabel(sprintf('Pre-IED %s change\nSpike-non spike difference',...
        pretty_name));
else
    
    ylabel(sprintf('Pre-IED %s\n%s change\nSpike-non spike difference',...
        metrics.time.freq(f).name,pretty_name));
end


set(gca,'fontsize',20)



%% Add p-value
sub_alpha = zeros(ntimes,1);
for tt = 1:last_before_rise-1
    num_left = ntimes-tt+1;
    num_sub_alpha = sum(pvals(tt:end) < 0.05/nfreq);
    if num_sub_alpha == num_left
        sub_alpha(tt) = 1;
    end
end

sig_times = find(sub_alpha == 1);
first_sig_time = sig_times(1);
fprintf('\nThe first significant AUC is %1.1f before the spike peak\n',...
    -times(first_sig_time));

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

endh = plot([mean_rise_spikes mean_rise_spikes],get(gca,'ylim'),'k--','linewidth',2);

legend(endh,'Visual rise','fontsize',20,'location','northwest')

if do_plot
print(gcf,[out_folder,sprintf('short_%s',met)],'-dpng');
end
    
    
end

