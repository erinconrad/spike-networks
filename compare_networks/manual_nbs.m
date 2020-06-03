function manual_nbs(whichPts,simple)

%% Parameters
plot_graphs = 0; % show graphs?
nperms = '5000';
nbs_method = 'Run NBS'; %'Run FDR' 'Run NBS'
NBS_test_threshold = '3.5';
test_method = 't-test';
alpha = '0.05';
graph_method = 'Extent'; %'Intensity' 'Extent'

% this specifically tests whether the group modeled by the first column is
% smaller than the group modeled by the second column (columns defined in
% design matrix). Given how I set up the design matrix, this tests whether
% the group in time i_time is more connected than the group in time 1
contrast = [-1 1]; % group 2 (later time) > group 1 (first time)

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

if simple == 1
    out_folder = [results_folder,'nbs_stats/simple/'];
elseif simple == 0
    out_folder = [results_folder,'nbs_stats/coherence/'];
end
if exist(out_folder,'dir') == 0
    mkdir(out_folder)
end

pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

sp_folder = [main_folder,'data/manual_spikes/'];
sp = load([sp_folder,'sp.mat']);
sp = sp.sp;

if isempty(whichPts) == 1
    whichPts = [];
    for i = 1:length(sp)
        if isempty(sp(i).name) == 0
            whichPts = [whichPts,i];
        end
    end
end



for whichPt = whichPts

    % Skip it if name is empty
    if isempty(pt(whichPt).name) == 1, continue; end
    name = pt(whichPt).name;
    
    if simple == 1
        adj_folder = [results_folder,'adj_mat/manual/adj_simple/'];
    elseif simple == 0
        adj_folder = [results_folder,'adj_mat/manual/adj_coherence/'];
    end

    %% Initialize structure to output stats
    nbs_stats.name = name;
    
    %% Load adjacency matrices 
    meta = load([adj_folder,name,'_adj.mat']);
    meta = meta.meta;
    
    % Get number of spikes, frequencies, times, channels
    nspikes = length(meta.spike);
    nfreq = length(meta.spike(1).adj);
    n_times = size(meta.spike(1).adj(1).adj,1);
    nchs = size(meta.spike(1).adj(1).adj,2);
    
    %% Get arrays of non-flattened adjacency matrices
    for i = 1:nfreq
        adj_avg(i).adj = zeros(n_times,nchs,nchs,nspikes);
        sim(i).p = zeros(n_times,1);
        sim(i).F = zeros(n_times,1);
        sim(i).indices = (1:n_times)';
        
         % First second is the control
        sim(i).p(1) = nan;
        sim(i).F(1) = nan;
    end
    
    % Loop over spikes
    for s = 1:length(meta.spike)
        
        for which_freq = 1:nfreq
        
            adj_all_t= meta.spike(s).adj(which_freq).adj;

            for t = 1:size(adj_all_t,1)
               
                % add to matrix
                adj_avg(which_freq).adj(t,:,:,s) = adj_all_t(t,:,:);

            end
            
        end
        
    end
    
    %% Remove columns with any nans or all zeros
    for which_freq = 1:nfreq
        adj_avg(which_freq).adj(:,:,:,any(isnan(adj_avg(which_freq).adj),[1 2 3])) = [];
    end
    
    for which_freq = 1:nfreq
        adj_avg(which_freq).adj(:,:,:,sum(adj_avg(which_freq).adj,[1 2 3])==0) = [];
    end
    
    %% Compare first time period to all additional time periods with permanova
    for which_freq = 1:nfreq        
        
        % Loop through subsequent times
        for i_time = 2:n_times
            
            fprintf('Doing time %d of %d, freq %d of %d\n',i_time,n_times,which_freq,nfreq);
            
            % Get design and adjacency matrices in the format needed for
            % NBS
            [mat,design] = transform_data_for_NBS(meta,which_freq,i_time);
                
            if 0
               % visualize the average adjacency matrices for the 2 conditions
               figure
               set(gcf,'position',[440 402 992 396])
               subplot(1,2,1)
               imagesc(nanmean(mat(:,:,1:size(mat,3)/2),3))
               title('Time 1')
               
               subplot(1,2,2)
               imagesc(nanmean(mat(:,:,size(mat,3)/2+1:end),3))
               title(sprintf('Time %d',i_time))
               pause
               close(gcf)
            end
            
            % Build a UI structure to feed into a command-line call to NBS
            UI.method.ui = nbs_method;
            UI.test.ui = test_method;
            UI.size.ui = graph_method;
            UI.thresh.ui = NBS_test_threshold;
            UI.perms.ui = nperms;
            UI.alpha.ui = alpha;
            UI.contrast.ui = contrast;
            UI.design.ui = design;
            UI.matrices.ui = mat;
            UI.exchange.ui=''; 
            UI.node_coor.ui = '';
            UI.node_label.ui = '';
            nbs = NBSrun(UI);
            
            % Visualize significant results
            fprintf('Number of significant sub-graphs for time %d: %d\n',...
                i_time,nbs.NBS.n);
            if plot_graphs == 1   
                for i = 1:nbs.NBS.n
                    figure
                    subplot(1,3,1)
                    % Plot average adjacency matrix for time 1
                    imagesc(squeeze(mean(adj_avg(which_freq).adj(1,:,:,:),4)))

                    subplot(1,3,2)
                    % Plot average adjacency matrix for time i_time
                    imagesc(squeeze(mean(adj_avg(which_freq).adj(i_time,:,:,:),4)))

                    subplot(1,3,3)
                    % Plot the connected components with significant changes
                    imagesc(nbs.NBS.con_mat{i})
                    pause
                    close(gcf)
                end 
            end
            
            % Fill up structure with stats
            nbs_stats.freq(which_freq).time(i_time).nbs = nbs.NBS;
            nbs_stats.freq(which_freq).time(i_time).parameters = nbs.UI;
            
        end
        
    end
        
%% Save the sim structure
save([out_folder,name,'_nbs.mat'],'nbs_stats');

end

end