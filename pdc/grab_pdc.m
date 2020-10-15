function stats = grab_pdc(time_window)


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

listing = dir([time_folder,'*.mat']);
stats.time.time_window = time_window;

all_names = {};
for l = 1:length(listing)
    fname = listing(l).name;
    pt_name = strsplit(fname,'_');
    pt_name = pt_name{1};
    [a,b] = ismember(pt_name,all_names);
    if a == 1
        pt_idx = b;
    else
        all_names = [all_names;pt_name];
        pt_idx = length(all_names);
    end

    
    meta = load([time_folder,fname]);
    meta = meta.meta;
    times = round(meta.spike(1).index_windows/meta.fs*1e2)/(1e2);
    
    for f = 1:length(meta.spike(1).adj)
        all_out_big = [];
        all_out_inv = [];
        for s = 1:length(meta.spike)
            big = meta.spike(s).biggest_dev;
            inv = meta.spike(s).is_sp_ch;

                adj = meta.spike(s).adj(f).adj;
                out_big = squeeze(sum(adj(:,big,:),3));
                out_inv = squeeze(mean(sum(adj(:,inv,:),3),2));
                
                all_out_big = [all_out_big,out_big];
                all_out_inv = [all_out_inv,out_inv];
        end
        
        stats.time.freq(f).pdc_big.pt(pt_idx).times = times;
        
        if ~contains(fname,'not')
            stats.time.freq(f).pdc_big.pt(pt_idx).spike.data = all_out_big';
            stats.time.freq(f).pdc_inv.pt(pt_idx).spike.data = all_out_inv';

        else
            stats.time.freq(f).pdc_big.pt(pt_idx).not.data = all_out_big';
            stats.time.freq(f).pdc_inv.pt(pt_idx).not.data = all_out_inv';
        end
    end
    
end

end