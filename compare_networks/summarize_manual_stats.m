function summarize_manual_stats(simple,time_window)

%{
This function only makes summary tables for the simple correlation case
%}

%% parameters
fs = 512;
alpha = 0.05;

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
bct_folder = locations.BCT;
addpath(genpath(bct_folder));
plot_folder = [results_folder,'plots/'];
time_text = sprintf('%1.1f/',time_window);

if simple == 1
    network_folder = [results_folder,'networks/manual/simple/',time_text];
    perm_folder = [results_folder,'perm_stats/simple/',time_text];
    nbs_folder = [results_folder,'nbs_stats/simple/',time_text];
    var_names_pt = {'Time','SignalDev','NBS','Perm','NS','GE'};
elseif simple == 0
    error('can only do for simple\n');
    network_folder = [results_folder,'networks/manual/coherence/'];
    perm_folder = [results_folder,'perm_stats/coherence/'];
    
end

% Load signal deviation file
sig_dev_folder = [results_folder,'signal_deviation/manual/'];
sig_dev = load([sig_dev_folder,'sig_dev.mat']);
sig_dev = sig_dev.sig_dev;


% Get the individual patient stat files
listing = dir([perm_folder,'*_perm.mat']);

% Initialize table of significant values
%sig_table = cell2table(cell(0,5),'VariableNames',{'Name','Freq','Time','F','p'});


% Loop through patients
pt_count = 0;
for i = 1:length(listing)
    
    
    % Load the permutation file
    name = listing(i).name;
    pt_name = strsplit(name,'_');
    pt_name = pt_name{1};
    sim = load([perm_folder,name]);
    sim = sim.sim;
    sim_p = sim.p;
    %sim_p_text = arrayfun(@(x) sprintf('%1.3f',x), sim_p,'UniformOutput',false);
    
    % Load the nbs stats file
    nbs_stats = load([nbs_folder,pt_name,'_nbs.mat']);
    nbs_stats = nbs_stats.nbs_stats;
    
    % Load network stats file
    net_stats = load([network_folder,pt_name,'_network_stats.mat']);
    metrics = net_stats.metrics;
    
    % Get an array of p values for nbs
    nbs_p = nan;
    for it = 2:length(nbs_stats.freq.time)
        nbs_curr = nbs_stats.freq.time(it).nbs;
        if nbs_curr.n == 0
            nbs_p = [nbs_p;nan];
        else
            % take the smallest p value of however many significant graphs
            % there are
            nbs_p = [nbs_p;min(nbs_curr.pval)];
        end
        
    end
    
    
    pt_count = pt_count + 1;
    
    nfreq = length(sim);
    if nfreq > 1, error('This function is only for simple correlation\n'); end
    
    % Get the appropriate times and signal deviation
    sig_dev_name = sig_dev(pt_count).name;
    if strcmp(sig_dev_name,pt_name) == 0
        pt_count = pt_count - 1;
        continue;
    end
    window = diff(sig_dev(pt_count).index_windows,1,2)/fs;
    window = window(1);
    times = 1:size(sig_dev(pt_count).index_windows,1);
    times = times*window-window;
    sig_dev_p = (sig_dev(pt_count).p);
    sig_dev_p_text = arrayfun(@(x) sprintf('%1.3f',x), sig_dev_p,'UniformOutput',false);
    
    
    
    % get metrics
    ns_p = metrics.metric(1).p;
    ge_p = metrics.metric(2).p;
    %ns_p_text = arrayfun(@(x) sprintf('%1.3f',x), ns_p,'UniformOutput',false);
    %ge_p_text = arrayfun(@(x) sprintf('%1.3f',x), ge_p,'UniformOutput',false);
    
    % Add asterixes when appropriate
    nbs_p_text = arrayfun(@(x,y) get_asterixes(x,y,alpha/(length(nbs_p)-1)),nbs_p,sig_dev_p','UniformOutput',false);
    sim_p_text = arrayfun(@(x,y) get_asterixes(x,y,alpha/(length(sim_p)-1)),sim_p,sig_dev_p','UniformOutput',false);
    ns_p_text = arrayfun(@(x,y) get_asterixes(x,y,alpha/(length(sim_p)-1)),ns_p',sig_dev_p','UniformOutput',false);
    ge_p_text=arrayfun(@(x,y) get_asterixes(x,y,alpha/(length(sim_p)-1)),ge_p',sig_dev_p','UniformOutput',false);
    
    % add stuff to table
    pt_name
    pt_table = table(times',sig_dev_p_text',nbs_p_text,sim_p_text,ns_p_text,...
        ge_p_text,'VariableNames',var_names_pt)
   
    
    
end


end


function out = get_asterixes(p,dev_p,alpha)
    %alpha = 0.05;
    out = sprintf('%1.3f',p);
    if p < alpha && dev_p > alpha
        out = [out,'***'];
    end

end