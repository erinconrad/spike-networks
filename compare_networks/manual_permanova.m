function manual_permanova(simple,time_window)

%{
This function compares adjacency matrices across time periods surrounding
the spikes. It flattens the matrices, and then does a multivariate
non-parametric test to compare all elements of the matrix across times.
%}

%% Parameters
nperms = 1e3;

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
    out_folder = [results_folder,'perm_stats/simple/',time_text];
    adj_folder = [results_folder,'adj_mat/manual/adj_simple/',time_text];
elseif simple == 0
    out_folder = [results_folder,'perm_stats/coherence/',time_text];
    adj_folder = [results_folder,'adj_mat/manual/adj_coherence/',time_text];
end

if exist(out_folder,'dir') == 0
    mkdir(out_folder)
end

pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

sp_folder = [main_folder,'data/manual_spikes/'];
sp = get_manual_times_from_excel;

%{
% old spikes
sp = load([sp_folder,'sp.mat']);
sp = sp.sp;
%}

listing = dir([adj_folder,'*_adj.mat']);



for j = 1:length(listing)

    filename = listing(j).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    
    
    %% Load adjacency matrices 
    meta = load([adj_folder,name,'_adj.mat']);
    meta = meta.meta;
    
    % Get number of spikes, frequencies, times, channels
    nspikes = length(meta.spike);
    nfreq = length(meta.spike(1).adj);
    n_times = size(meta.spike(1).adj(1).adj,1);
    nchs = size(meta.spike(1).adj(1).adj,2);
    
    %% Get arrays of flattened adjacency matrices
    for i = 1:nfreq

        adj_avg(i).adj = nan(n_times,nchs*(nchs-1)/2,nspikes);
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
                % flatten matrix
                adj_out = flatten_or_expand_adj(squeeze(adj_all_t(t,:,:)));

                % add to matrix
                adj_avg(which_freq).adj(t,:,s) = adj_out;

            end
            
        end
        
    end
    
    %% Remove columns with any nans or all zeros
    for which_freq = 1:nfreq
        adj_avg(which_freq).adj(:,:,any(isnan(adj_avg(which_freq).adj),[1 2])) = [];
    end

    for which_freq = 1:nfreq
        adj_avg(which_freq).adj(:,:,sum(adj_avg(which_freq).adj,3)==0) = [];
    end
    
    
    %% Compare first time period to all additional time periods with permanova
    for which_freq = 1:nfreq
        % Get the first second (this will be an nch*(nch-1)/2 x nspike
        % matrix)
        first_sec = squeeze(adj_avg(which_freq).adj(1,:,:));

        % reshape it for permanova (so now it's nspike x nch*(nch-1)/2)
        first_sec = first_sec';

        % make a nspike x 1 array with time 1
        first_sec_time = ones(size(first_sec,1),1);
        
        % Loop through subsequent times
        for i_time = 2:n_times
            
            fprintf('Doing time %d of %d, freq %d of %d\n',i_time,n_times,which_freq,nfreq);
            
            % Get current time
            test_sec = squeeze(adj_avg(which_freq).adj(i_time,:,:));
                
            % Reshape it for perMANOVA
            test_sec = test_sec';
            next_sec_time = 2*ones(size(test_sec,1),1);
            
            % Do perMANOVA
            
            % Generate the dissimilarity matrix. When I supply this without
            % arguments it defaults to euclidean distance, which is to say
            % the sum of the square difference between each feature of the
            % nspikex1 vector between time period 1 and 2. This generates
            % an nspikes*2 x nspikes*2 dissimilarity matrix, basically how
            % different each spike is from each other spike. If first_sec
            % and test_sec are very different, then I would expect the RUQ
            % and LLQ of the dis matrix to be very bright relative to the
            % LUQ and RLQ (saying that the spikes in different time periods
            % are more different from each other than they are from spikes
            % in the same time period).
            dis = f_dis([first_sec;test_sec]);
            
            % Get an F statistic and p-value on the dissimilarity between
            % time period 1 and 2 using permutation test on the
            % dissimilarity matrix
            result = f_permanova(dis,[first_sec_time;next_sec_time],nperms);
            
            sim(which_freq).p(i_time) = result.p;
            sim(which_freq).F(i_time) = result.F;
            
            
            if 0
                % Plot the dissimilarity matrix. Each i,j element shows how
                % dissimilar time i is from time j (brighter colors mean
                % more dissimilar). Remember that the first half of the
                % elements are from time period 1 and the second half are
                % from time period 2.
                figure
                imagesc(dis)
                xticks(1:size(dis,1))
                yticks(1:size(dis,1))
                title(sprintf('Time period 1 vs %d, p = %1.3f',i_time,result.p))
                pause
                close(gcf)
            end
            
        end
        
    end
        
%% Save the sim structure
save([out_folder,name,'_perm.mat'],'sim');

end

end