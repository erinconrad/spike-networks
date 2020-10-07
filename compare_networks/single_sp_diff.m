function single_sp_diff(overwrite,simple,time_window,not_a_spike)

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
time_text = sprintf('%1.1f/',time_window);

if simple == 1
    out_folder = [results_folder,'net_diff_stats/simple/',time_text];
    adj_folder = [results_folder,'adj_mat/manual/adj_simple/',time_text];
elseif simple == 0
    out_folder = [results_folder,'net_diff_stats/coherence/',time_text];
    adj_folder = [results_folder,'adj_mat/manual/adj_coherence/',time_text];
end

if not_a_spike
    not_a_spike_text = '_not_spike_';
else
    not_a_spike_text = '_';
end

if exist(out_folder,'dir') == 0
    mkdir(out_folder)
end

pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

listing = dir([adj_folder,'*_adj.mat']);


finished_names = {};
for j = 1:length(listing)

    filename = listing(j).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    
    if overwrite == 0
        if exist([out_folder,name,not_a_spike_text,'perm.mat'],'file') ~= 0
            fprintf('%s already exists, skipping...\n',name);
            continue;
        end
    end
    
    if ismember(name,finished_names)
        continue
    end
    
    fprintf('\nDoing %s\n',name);
    
    %% Load adjacency matrices 
    meta = load([adj_folder,name,not_a_spike_text,'adj.mat']);
    meta = meta.meta;
    
    % Get number of spikes, frequencies, times, channels
    nspikes = length(meta.spike);
    nfreq = length(meta.spike(1).adj);
    n_times = size(meta.spike(1).adj(1).adj,1);
    nchs = size(meta.spike(1).adj(1).adj,2);
    
    %% Get arrays of flattened adjacency matrices
    for i = 1:nfreq
        sim(i).pt_name = name;
        adj_avg(i).adj = nan(n_times,nchs*(nchs-1)/2,nspikes);
        if isfield(meta.spike(1).adj(i),'name') == 1
            sim(i).name = meta.spike(1).adj(i).name;
        else
            sim(i).name = 'correlation';
        end
        sim(i).F = zeros(n_times,nspikes);
        sim(i).indices = (1:n_times)';
        
        % First second is the control
        sim(i).F(1,:) = nan(1,nspikes);
        
        
        sim(i).index_windows = meta.spike(1).index_windows;
        sim(i).fs = meta.fs;
        sim(i).times = round(sim(i).index_windows(:,1)/sim(i).fs*1e2)/(1e2);
       
    end
    
    % Loop over spikes
    for s = 1:length(meta.spike)
        
        for which_freq = 1:nfreq
        
            adj_all_t= meta.spike(s).adj(which_freq).adj;

            for t = 1:size(adj_all_t,1)
                % flatten matrix
                adj_out = flatten_or_expand_adj(squeeze(adj_all_t(t,:,:)));

                % add to matrix
                adj_avg(which_freq).adj(t,:,s) = adj_out;

            end
            
            adj_1 = squeeze(adj_avg(which_freq).adj(1,:,s));
            % compare first time to all subsequent times
            for t = 2:size(adj_all_t,1)
                adj_t = squeeze(adj_avg(which_freq).adj(t,:,s));
                net_diff = sum((adj_1-adj_t).^2); %sum square diff
                sim(which_freq).F(t,s) = net_diff;
            end
            
        end
        
    end
    
    % Optional plot
    if 0
        for f = 1:nfreq
            figure
            times = (1:size(sim(f).F,1))';
            plot(sim(f).F,'o')
            hold on
            for tt = 1:length(times)
                plot([times(tt)-0.2 times(tt)+0.2],...
                    [mean(sim(f).F(tt,:)) mean(sim(f).F(tt,:))])
            end
            xlabel('Time window')
            ylabel('Net diff')
            title(sprintf('%s %s',name,sim(i).name))
            pause
            close(gcf)
        end
    end
    
    %% PCA
    
    for f = 1:length(adj_avg)
        old_adj = adj_avg(f).adj;
        new_adj = zeros(size(old_adj,1)*size(old_adj,3),size(old_adj,2));
        time_idx = zeros(size(old_adj,1)*size(old_adj,3),1);
        
        % Concatenate different times together
        for t = 1:size(old_adj,1)
            new_adj(size(old_adj,3)*(t-1)+1:size(old_adj,3)*t,:) = ...
                (squeeze(old_adj(t,:,:)))';
            time_idx(size(old_adj,3)*(t-1)+1:size(old_adj,3)*t) = ...
                ones(size(old_adj,3),1)*t;
        end
        
        % Do pca
        [coeff,score,latent] = pca(new_adj);
        
        % Plot the scores over time periods for the first three components
        if 0
        for c = 1:3
            
            figure
            for t = 1:size(old_adj,1)
                % get the indices of the networks in those time windows
                curr_idx = time_idx == t;
                
                % get the scores of this component for those time windows
                curr_scores = score(curr_idx,c);
                
                plot(t,curr_scores,'ko');
                hold on
                plot([t-0.2 t+0.2],[mean(curr_scores) mean(curr_scores)],'k-');
            end
        end
        end
        
        % Save the scores of the first component
        sim(f).score = score(:,1);
        sim(f).time_idx = time_idx;
    end
        
    %% Save the sim structure
    save([out_folder,name,not_a_spike_text,'perm.mat'],'sim');

    %% Add to finished names
    finished_names = [finished_names,name];
end

end