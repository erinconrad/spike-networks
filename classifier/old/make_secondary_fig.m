function make_secondary_fig(stats,windows,met,which_pre_rise,which_freq)

if strcmp(met,'sd')
    met_text = 'power';
elseif strcmp(met,'ers')
    met_text = [newline,'frequency-specific power'];
elseif strcmp(met,'F')
    met_text = 'network difference';
else
    met_text = met;
end

locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
out_folder = [results_folder,'plots/'];

if length(windows) > 1, error('Only one time window'); end

network_count = length(stats);
n_freq_abs = length(stats(1).time(1).freq);

%% Initialize figure
figure
set(gcf,'position',[1 100 800 500])


% Find the right window
for t = 1:length(stats(1).time)
    if ismember(stats(1).time(t).time_window,windows), break; end
end


%% Plot the thing
n = 1;

net_name = stats(n).name;
tcount = 1;

times = stats(n).time(t).freq(1).(met).pt(1).times;
nfreq = length(stats(n).time(t).freq);

f = which_freq;

% Get appropriate subplot
if strcmp(net_name,'coherence') == 1
    column_add = 1;
else
    column_add = 0;
end
sp = (n_freq_abs+1)*(tcount) + f + column_add;

dat_sp = stats(n).time(t).freq(f).(met).all_z_spike;
dat_not = stats(n).time(t).freq(f).(met).all_z_not;

% Plot them
z_range = [0,0];
[z_range,prettyp] = specific_plot(dat_sp,dat_not,0,z_range,1,times,n_freq_abs+1);

xlabel('Time relative to spike peak (s)')
ylabel(sprintf('Normalized %s',met_text))    
xl = get(gca,'xlim');
yl = get(gca,'ylim');

text((xl(1)+xl(2))/2,yl(1)+0.9*(yl(2)-yl(1)),...
    sprintf('%s',prettyp),'fontsize',20,...
    'horizontalalignment','center')




print(gcf,[out_folder,'isolated_',sprintf('%s_%d',met,which_pre_rise)],'-depsc');

end