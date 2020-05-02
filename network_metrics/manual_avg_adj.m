function manual_avg_adj(whichPts,simple)

%% Parameters
which_times = [8 9 10 11 12 13 14];

% 1 = alpha/theta; 2 = beta, 3 = low gamma, 4 = high gamma, 5 = ultra high, 6 = broadband
freq_text = {'alpha/theta','beta','low\ngamma','high\ngamma','ultra high\ngamma','broadband'};


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
if exist(plot_folder,'dir') == 0
    mkdir(plot_folder);
end


pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;


for whichPt = whichPts
    
    %% Prep patient
    
    % Skip it if name is empty
    if isempty(pt(whichPt).name) == 1, continue; end
    name = pt(whichPt).name;
    
    if simple == 1
        adj_folder = [results_folder,'adj_mat/manual/adj_simple/'];
    elseif simple == 0
        adj_folder = [results_folder,'adj_mat/manual/adj_coherence/'];
    end

    
    %% Load adjacency matrices and calculate metrics
    meta = load([adj_folder,name,'_adj.mat']);
    meta = meta.meta;
    
    %% Prep avg adjacency matrices
    % Load the first adjacency matrix to get the size in order to
    % initialize one of the arrays
    nfreq = length(meta.spike(1).adj);
      
    for i = 1:nfreq
        adj_avg(i).adj = zeros(size(meta.spike(1).adj(1).adj));
    end
    
    % Initialize spike count
    s_count = 0;

    % Loop through spikes
    for s = 1:length(meta.spike)

        for which_freq = 1:nfreq
            
            if sum(sum(sum(isnan(meta.spike(s).adj(which_freq).adj)))) > 0
                continue
            end
            adj_all_t= meta.spike(s).adj(which_freq).adj; 
            
            adj_avg(which_freq).adj = adj_avg(which_freq).adj + adj_all_t;
        end

        s_count = s_count + 1;

    end

end
    
    %% Divide by number of spikes to get average and calculate global metrics

    for which_freq = 1:length(adj_avg)
        adj_avg(which_freq).adj = adj_avg(which_freq).adj/s_count;
        
        if sum(sum(sum(isnan(adj_avg(which_freq).adj)))) > 0
            continue
        end
        

    end
    
   
    %% Plot the average adjacency matrices
    figure
    set(gcf,'position',[0 0 (length(which_times)-1)*150+300 (nfreq-1)*200+200])
    [ha, pos] = tight_subplot(nfreq,length(which_times), [0.07 0.02], [0.12 0.15], [0.05 0.01]);
    count = 0;
    for f = 1:nfreq
        for it = 1:length(which_times)
            count = count + 1;
            axes(ha(count))
            
            imagesc(squeeze(adj_avg(f).adj(which_times(it),:,:)))
            title(sprintf('Time %d',which_times(it)))
            xticklabels([])
            yticklabels([])
            
            if it == 1
                if simple == 1
                    ylabel('Pairwise correlation')
                elseif simple == 0
                    ylabel(sprintf('%s',freq_text{f}))
                end
            end
            
            set(gca,'fontsize',20)
        end
    end
    
    

end