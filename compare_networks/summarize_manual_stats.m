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
    out_folder = [results_folder,'stats_sum/simple/',time_text];
elseif simple == 0
    network_folder = [results_folder,'networks/manual/coherence/',time_text];
    perm_folder = [results_folder,'perm_stats/coherence/',time_text];
    nbs_folder = [results_folder,'nbs_stats/coherence/',time_text];
    adj_folder = [results_folder,'adj_mat/manual/adj_coherence/',time_text];
    out_folder = [results_folder,'stats_sum/coherence/',time_text];
end
var_names_pt = {'Time','SignalDev','NBS','Perm','NS','GE'};

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

% Load signal deviation file
sig_dev_folder = [results_folder,'signal_deviation/manual/',time_text];
sig_dev = load([sig_dev_folder,'sig_dev.mat']);
sig_dev = sig_dev.sig_dev;


% Get the individual patient stat files
listing = dir([perm_folder,'*_perm.mat']);

% Initialize table of significant values
%sig_table = cell2table(cell(0,5),'VariableNames',{'Name','Freq','Time','F','p'});

% Loop through patients
for i = 1:length(listing)
    
    % Get pt name
    name = listing(i).name;
    pt_name = strsplit(name,'_');
    pt_name = pt_name{1};
    
    % Load the adjacency matrix to get frequency bands (for coherence
    % measurements)
    if simple == 0
        meta = load([adj_folder,pt_name,'_adj.mat']);
        meta = meta.meta;
        freq_text = {};
        for f = 1:length(meta.spike(1).adj)
            freq_text = [freq_text,meta.spike(1).adj(f).name];
        end
        nf = length(meta.spike(1).adj);
        nt = size(meta.spike(1).index_windows,1);
    else
        nf = 1;
        
    end
    
    
    % Load the permutation file and get perm stats
    sim = load([perm_folder,name]);
    sim = sim.sim;
    if simple == 1
        sim_p = sim.p;
        nt = length(sim_p);
    else
        sim_p = nan(nt,nf);
        for f = 1:length(sim)
            sim_p(:,f) = sim(f).p;
        end
    end
         
    
    % Load the nbs stats file
    nbs_stats = load([nbs_folder,pt_name,'_nbs.mat']);
    nbs_stats = nbs_stats.nbs_stats;
    
    % Get an array of p values for nbs
    if simple == 1
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
    else
        nbs_p = nan(nt,nf);
        for f = 1:nf
            for t = 2:nt
                nbs_curr = nbs_stats.freq(f).time(t).nbs;
                if nbs_curr.n == 0
                    nbs_p(t,f) = nan;
                else
                    nbs_p(t,f) = min(nbs_curr.pval);
                end
            end
        end
    end
     
    
    % Find the appropriate pt in signal deviation
    sig_dev_pt = nan;
    for p = 1:length(sig_dev)
        if strcmp(sig_dev(p).name,pt_name) == 1
            sig_dev_pt = p;
            break
        end
    end
    if isnan(sig_dev_pt) == 1, error('cannot find pt\n'); end
    
    sig_dev_name = sig_dev(sig_dev_pt).name;
    % Get the appropriate times and signal deviation
    window = diff(sig_dev(sig_dev_pt).index_windows,1,2)/fs;
    window = window(1);
    times = 1:size(sig_dev(sig_dev_pt).index_windows,1);
    times = times*window-window;
    sig_dev_p = (sig_dev(sig_dev_pt).p);
    sig_dev_p_text = arrayfun(@(x) sprintf('%1.3f',x), sig_dev_p,'UniformOutput',false);
    
    
    % Load network stats file
    if exist([network_folder,pt_name,'_network_stats.mat'],'file') == 0
        net_stats = [];
    else
        net_stats = load([network_folder,pt_name,'_network_stats.mat']);
    end
    
    if simple == 1
        n_comparisons = nt-1;
    else
        n_comparisons = (nt-1)*nf;
    end
    
    % get metrics
    if isempty(net_stats) == 0
        metrics = net_stats.metrics;
        if simple == 1
            ns_p = metrics.metric(1).p;
            ge_p = metrics.metric(2).p;
            ns_p(1) = nan; % first time currently set to 0 and not nan
            ge_p(1) = nan;
        else
            ns_p = metrics.metric(1).p;
            ge_p = metrics.metric(2).p;
            ns_p(:,1) = nan; % first time currently set to 0 and not nan
            ge_p(:,1) = nan;
        end
    end
    
    
    % Add asterixes when appropriate
    if isempty(net_stats) == 0
    ns_p_text = arrayfun(@(x,y) get_asterixes(x,y,alpha/n_comparisons),ns_p',...
        repmat(sig_dev_p',1,nf),'UniformOutput',false);
    ge_p_text=arrayfun(@(x,y) get_asterixes(x,y,alpha/n_comparisons),ge_p',...
        repmat(sig_dev_p',1,nf),'UniformOutput',false);
    end
    nbs_p_text = arrayfun(@(x,y) get_asterixes(x,y,alpha/n_comparisons),nbs_p,...
        repmat(sig_dev_p',1,nf),'UniformOutput',false);
    sim_p_text = arrayfun(@(x,y) get_asterixes(x,y,alpha/n_comparisons),sim_p,...
        repmat(sig_dev_p',1,nf),'UniformOutput',false);

    
    
    % add stuff to table
    pt_name
    if isempty(net_stats) == 0
        if simple == 1
            pt_table = table(times',sig_dev_p_text',nbs_p_text,sim_p_text,ns_p_text,...
                ge_p_text,'VariableNames',var_names_pt)
            all_tables(i).name = pt_name;
            all_tables(i).freq.name = 'na';
            all_tables(i).freq.table = pt_table;
        else
            for f = 1:nf
                freq_text{f}
                pt_table = table(times',sig_dev_p_text',nbs_p_text(:,f),...
                    sim_p_text(:,f),ns_p_text(:,f),...
                ge_p_text(:,f),'VariableNames',var_names_pt)
                all_tables(i).name = pt_name;
                all_tables(i).freq(f).name = freq_text{f};
                all_tables(i).freq(f).table = pt_table;
            end
        end
    end
   
    
    
end


save([out_folder,'stats_sum.mat'],'all_tables');

end


function out = get_asterixes(p,dev_p,alpha)
    %alpha = 0.05;
    out = sprintf('%1.3f',p);
    if p < alpha && dev_p > alpha
        out = [out,'***'];
    end

end