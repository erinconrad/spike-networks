function [percs,all,total_soz_n,all_main,all_soz_chs] = soz_info

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
eeg_folder = [main_folder,'results/eeg_data/'];
data_folder = [main_folder,'data/'];

bct_folder = locations.BCT;
addpath(genpath(bct_folder));
out_folder = [results_folder,'plots/'];
bad_folder = [results_folder,'bad_spikes/'];
seq_folder = [results_folder,'seq_data/'];

pt = load([data_folder,'spike_structures/pt.mat']);
pt = pt.pt;

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

pt_listing = dir([eeg_folder,'*.mat']);

all_big_in_soz = [];
all_lead_in_soz = [];

perc_big_in_soz = [];
perc_lead_in_soz = [];

total_soz_n = [];
all_main = [];
all_soz_chs = {};

        % loop through pts
for i = 1:length(pt_listing)
    
    fname = pt_listing(i).name;
    pt_name = strsplit(fname,'_');
    pt_name = pt_name{1};
    
    spike = load([eeg_folder,pt_name,'_eeg.mat']);
    spike = spike.spike;
    
    % Load bad file
    bad = load([bad_folder,pt_name,'_rise.mat']);
    bad = bad.bad;
    
    % Load seq file
    seq = load([seq_folder,pt_name,'_seq.mat']);
    seq = seq.seq;
    
    % Get soz
    soz_chs = get_soz_chs(pt,pt_name);
    
    % identities of soz channels
    all_soz_chs = [all_soz_chs;soz_chs];
    
    % number of soz channels and total chs
    total_soz_n = [total_soz_n;length(soz_chs),length(spike(1).chLabels)];
    
    out(i).name = pt_name;
    out(i).in_soz = [];
    out(i).all = [];
    
    for s = 1:length(spike)
        
        % Skip if bad
        if bad.spike(s).bad_spike == 1
            continue
        end
        
        % Biggest dev
        big = spike(s).biggest_dev;
        big_in_soz = ismember(big,soz_chs);
        
        % involved
        inv = find(spike(s).involved);
        inv_in_soz = ~isempty(intersect(inv,soz_chs));
        
        % seq
        if isempty(seq(s).seq)
            sq = [];
            lead = nan;
            lead_in_soz = nan;
            
        else
            sq = seq(s).seq(:,1);
            lead = seq(s).first_ch;
            lead_in_soz = ismember(lead,soz_chs);
        end
        
        out(i).in_soz = [out(i).in_soz;big_in_soz, inv_in_soz, lead_in_soz];
        out(i).all = [out(i).all;big, lead];
        
    end
    
    main_big = mode(out(i).all(:,1));
    main_lead = mode(out(i).all(:,2));
    all_main = [all_main;main_lead main_big];
    
    perc_big_in_soz = [perc_big_in_soz;sum(out(i).in_soz(:,1))/length(out(i).in_soz(:,1))];
    perc_lead_in_soz = [perc_lead_in_soz;nansum(out(i).in_soz(:,3))/sum(~isnan(out(i).in_soz(:,3)))];
    
    
    main_big_in_soz = ismember(main_big,soz_chs);
    main_lead_in_soz = ismember(main_lead,soz_chs);
    
    out(i).main_big_in_soz = main_big_in_soz;
    out(i).main_lead_in_soz = main_lead_in_soz;
    
    all_big_in_soz = [all_big_in_soz;main_big_in_soz];
    all_lead_in_soz = [all_lead_in_soz;main_lead_in_soz];
end

percs = [perc_lead_in_soz,perc_big_in_soz];
all = [all_lead_in_soz,all_big_in_soz];
end