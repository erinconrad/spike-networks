function plot_pdc

time_window = 0.5;
met = 'pdc_big';

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
eeg_folder = [main_folder,'results/eeg_data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));

bct_folder = locations.BCT;
addpath(genpath(bct_folder));
out_folder = [results_folder,'plots/'];
sig_dev_folder = [results_folder,'signal_deviation/manual/'];
perm_folder = [results_folder,'perm_stats/'];
ers_folder = [results_folder,'ers/'];
ns_folder = [results_folder,'metrics/manual/'];
sp_diff_folder = [results_folder,'net_diff_stats/'];
pdc_folder = [results_folder,'pdc/'];

time_folder = [pdc_folder,sprintf('%1.1f/',time_window)];

%% grab pdc
stats = load([time_folder,'all_pdc.mat']);
stats = stats.stats;
nfreq = length(stats.time.freq);



for f = 1:nfreq
    curr_met = stats.time.freq(f).(met);
    for sp = {'spike','not'}
        curr_met.(sp{1}).all_z = [];
        for p = 1:length(curr_met.pt)

            all_z = [];
            for s = 1:size(curr_met.pt(p).(sp{1}).data,1)
                data = curr_met.pt(p).(sp{1}).data(s,:);
                z = (data-nanmean(data))/nanstd(data);
                all_z = [all_z;z];
            end
            mean_z = nanmean(all_z,1);
            curr_met.(sp{1}).all_z = [curr_met.(sp{1}).all_z;mean_z];
            stats.time.freq(f).(met) = curr_met;
        end
    end
        
end

set(gcf,'position',[1 100 1500 300])
[ha, ~] = tight_subplot(1, nfreq, [0.08 0.01], [0.2 0.12], [0.07 0.005]);
for f = 1:nfreq
    axes(ha(f))
    times = stats.time.freq(f).(met).pt(1).times(:,1);
    dat_sp = stats.time.freq(f).(met).spike.all_z;
    dat_not = stats.time.freq(f).(met).not.all_z;
    mean_sp = nanmean(dat_sp,1);
    mean_not = nanmean(dat_not,1);
    
    ste_sp = nanstd(dat_sp,0,1)./sqrt(sum(~isnan(dat_sp),1));
    ste_not = nanstd(dat_not,0,1)./sqrt(sum(~isnan(dat_not),1));
    errorbar(times,mean_sp,ste_sp)
    hold on
    errorbar(times,mean_not,ste_not)
end


end

