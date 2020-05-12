function dev_vs_network(whichPts,small)

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));

spike_times_file = [data_folder,'spike_times/times.mat'];
pt_file = [data_folder,'spike_structures/pt.mat'];
pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;
times = load(spike_times_file); % will result in a structure called "out"
times = times.out;

alpha = 0.05/21; %Bonferroni correction for 21 comparisons.
spike_periods = [11 12];


freq_names = {'alpha_theta','beta','low_gamma',...
    'high_gamma','ultra_high','broadband'};

for whichPt = whichPts
    
    % Skip it if name is empty
    if isempty(times(whichPt).name) == 1, continue; end
    
    
    name = times(whichPt).name;

    dev_folder = [results_folder,'sig_dev/',name,'/'];
    pt_folder = [results_folder,name,'/'];
    
    if exist(dev_folder,'dir') == 0, continue; end
    
    if small == 3
        net_folder = [pt_folder,'stats_simple/'];
    elseif small == 4
        net_folder = [pt_folder,'stats_coherence/'];
    end
    
    if exist(net_folder,'dir') == 0, continue; end
    
    % Load signal deviation significance
    dev = load([dev_folder,'sig_dev.mat']);
    dev_p = dev.sig_dev.p;
    
    % Load network significance
    net = load([net_folder,'sim_permanova.mat']);
    
    % Find significant signal deviations
    dev_sig = find(dev_p < alpha);
    
   
    
    fprintf('\n');
 
    % NOT BONFERRONI CORRECTING FOR MULTIPLE FREQ BANDS
    for which_freq = 1:length(net.sim)
        
        net_p = net.sim(which_freq).p;
        
        if small == 3
            fprintf('%s deviation p and network p:\n',name);
        elseif small == 4
            fprintf('%s deviation p and network p, %s:\n',name,freq_names{which_freq});
        end
        [(1:22)',dev_p',net_p]

        %% Say if the spike time periods (11 and 12 are significant)
        net_sig = find(net_p < alpha);
        sp_period_sig = ismember(spike_periods,net_sig);
        sp_period_sig = spike_periods(sp_period_sig);
        if isempty(sp_period_sig) == 1
            if small == 3
                fprintf('%s: NO SIGNIFICANT SPIKE NETWORK CHANGE.\n',name);
            elseif small == 4
                fprintf('%s %s: NO SIGNIFICANT SPIKE NETWORK CHANGE.\n',name,freq_names{which_freq});
            end
        else
            for i = 1:length(sp_period_sig)
                if small == 3
                    fprintf('%s: significant spike network change at period %d.\n',name,sp_period_sig(i));
                elseif small == 4
                    fprintf('%s %s: significant spike network change at period %d.\n',...
                        name,freq_names{which_freq},sp_period_sig(i));
                end
            end
        end

        %% Say if pre- or post-spike time periods are significant
        
        % Remove spike periods
        net_sig(ismember(net_sig,spike_periods)) = [];
        
        % Find those significant for network change but NOT signal
        % deviation
        sig_only_net = ~ismember(net_sig,dev_sig);
        sig_only_net = net_sig(sig_only_net);
        
        if isempty(sig_only_net) == 1
            if small == 3
                fprintf('%s: no significant pre- or post-spike network change.\n',name);
            elseif small == 4
                fprintf('%s %s: no significant pre- or post-spike network change.\n',name, freq_names{which_freq});
            end
        else
            for i = 1:length(sig_only_net)
                if small == 3
                    fprintf('%s: significant peri-spike network change for %d.\n',name,sig_only_net(i));
                elseif small == 4
                    fprintf('%s %s: significant peri-spike network change for %d.\n',...
                        name,freq_names{which_freq},sig_only_net(i));
                end
            end
        end
        
    end
    
    

end

end