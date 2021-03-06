function nbs_change_fig(graph_method)

%{
I am not sure how to plot the NBS statistics, not sure what to show other
than a p-value which is kind of lame. Thinking I will just show the
permanova statistic.
%}

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
eeg_folder = [main_folder,'results/eeg_data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
bct_folder = locations.BCT;
addpath(genpath(bct_folder));
out_folder = [results_folder,'plots/'];
sig_dev_folder = [results_folder,'signal_deviation/manual/'];

freq_names = {'delta','theta','alpha','beta','low_gamma',...
    'high_gamma','ultra_high','broadband'};

if strcmp(graph_method,'Extent') == 1
    main_folder_text = 'nbs_stats';
elseif strcmp(graph_method,'Intensity') == 1
    main_folder_text = 'nbs_stats_intensity';
end

perm_folder = [results_folder,main_folder_text,'/'];


if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

% Load spike file for one patient to get the surround time
spike = load([eeg_folder,'HUP074_eeg.mat']);
spike = spike.spike;
surround_time = spike(1).surround_time;


% Loop through network types
listing = dir(perm_folder);
network_count = 0;
n_freq_abs = 0;
max_F = 0;
max_z = 0;
min_z = 0;
network_names = {};
for l = 1:length(listing)
    name= listing(l).name;
    
    % Skip if . or ..
    if strcmp(name,'.') == 1 || strcmp(name,'..') == 1 || strcmp(name,'.DS_Store') == 1
        continue
    end
    network_names = [network_names;name];
    
    % Skip if not a directory
    if listing(l).isdir == 0, continue; end
    
    network_count = network_count + 1;
    stats(network_count).name = name;
    
    network_folder = [perm_folder,name,'/'];
    
    % Loop through time scales
    time_listing = dir(network_folder);
    time_count = 0;
    time_names = {};
    
    for k = 1:length(time_listing)
        time_name= time_listing(k).name;
        time_window = str2num(time_name);
        
        % Skip if . or ..
        if strcmp(time_name,'.') == 1 || strcmp(time_name,'..') == 1 || strcmp(time_name,'.DS_Store') == 1
            continue
        end
        time_names = [time_names;time_name];

        % Skip if not a directory
        if time_listing(k).isdir == 0, continue; end
        
        time_count = time_count + 1;
        stats(network_count).time(time_count).name = time_name;
        stats(network_count).time(time_count).time_window = time_window;
        time_folder = [network_folder,time_name,'/'];
        
        pt_listing = dir([time_folder,'*.mat']);
        
        % load one to get nfreq
        sim = load([time_folder,pt_listing(1).name]);
        sim = sim.nbs_stats;
        nfreq = length(sim.freq);
        if n_freq_abs < nfreq
            n_freq_abs = nfreq;
        end
        
       
        if isfield(sim,'index_windows')
            
            stats(network_count).time(time_count).index_windows = sim.index_windows;
            stats(network_count).time(time_count).fs = sim.fs;
            for f = 1:nfreq
                stats(network_count).time(time_count).freq(f).n_all = ...
                    nan(length(pt_listing),size(sim.index_windows,1));

                stats(network_count).time(time_count).freq(f).p_all = ...
                    ones(length(pt_listing),size(sim.index_windows,1));
            end
        else
            for f = 1:nfreq
                stats(network_count).time(time_count).freq(f).n_all = ...
                    nan(length(pt_listing),surround_time*2/time_window);

                stats(network_count).time(time_count).freq(f).p_all = ...
                    ones(length(pt_listing),surround_time*2/time_window);

            end
        end
        
        pt_names = {};
        % loop through pts
        for i = 1:length(pt_listing)
            
            pt_name = pt_listing(i).name;
            pt_name_pt = strsplit(pt_name,'_');
            pt_name_pt = pt_name_pt{1};
            pt_names = [pt_names;pt_name_pt];
            
            % load pt file
            sim = load([time_folder,pt_name]);
            sim = sim.nbs_stats;
            
            
            for f = 1:nfreq
                
                
                if isfield(sim.freq(f),'name')
                    stats(network_count).time(time_count).freq(f).name = sim.freq(f).name;
                else
                    fprintf('\nWarning, frequency not specified, using backup info.\n');
                    if nfreq == 1
                        stats(network_count).time(time_count).freq(f).name = 'correlation';
                    else
                        stats(network_count).time(time_count).freq(f).name = freq_names{f};
                    end
                end
                
                    for tt = 2:length(sim.freq(f).time)
                        stats(network_count).time(time_count).freq(f).n_all(i,tt) = sim.freq(f).time(tt).nbs.n;
                        if isempty(sim.freq(f).time(tt).nbs.pval) == 1
                            stats(network_count).time(time_count).freq(f).p_all(i,tt) = 1;
                        else
                            temp_p = min(sim.freq(f).time(tt).nbs.pval);
                            if temp_p == 0
                                temp_p = 1/str2double(sim.freq(f).time(tt).parameters.nperms);
                            end
                            stats(network_count).time(time_count).freq(f).p_all(i,tt) = temp_p;
                        end
                    
                    end
            
            end
            
            
        end
        
    end
 
end

%% Plot
% Initialize figure
%{
nfreq (coherence) + 1 (simple) columns and 2 (2 time scales) rows
%}
figure
set(gcf,'position',[1 100 1399 500])
[ha, pos] = tight_subplot(time_count, n_freq_abs+1, [0.01 0.01], [0.1 0.05], [0.05 0.01]);

for n = 1:network_count

    net_name = stats(n).name;
    
    for t = 1:time_count
        
        if t>size(stats(n).time), continue; end
        
        % change times for x axis
        nchunks = size(stats(n).time(t).freq(1).n_all,2);
        
        if isfield(stats(n).time(t),'index_windows') && isempty(stats(n).time(t).index_windows) == 0
            temp_times = stats(n).time(t).index_windows(:,1)/stats(n).time(t).fs-3;
            times = realign_times(temp_times,surround_time);
        else

            times = realign_times(nchunks,surround_time);
        end
        
        
        nfreq = length(stats(n).time(t).freq);
        for f = 1:nfreq
            
            % Get appropriate subplot
            if strcmp(net_name,'coherence') == 1
                column_add = 1;
            else
                column_add = 0;
            end
            % this adds the number of frequencies + 1 if it's on the 2nd
            % time point (to move down a row), and it adds which frequency
            % (which is 1 if simple) and adds 1 if coherence, to start with
            % the 2nd column for coherence
            sp = (n_freq_abs+1)*(t-1) + f + column_add;
            axes(ha(sp));
          
            n_curr = stats(n).time(t).freq(f).n_all;

            % plot mean F across patients
            for tt = 1:size(n_curr,2)
                
                curr_p_vals = stats(n).time(t).freq(f).p_all(:,tt);
                comb_p = fisher_p_value(curr_p_vals);
                text_out = get_asterisks(comb_p,(nchunks-1)*(n_freq_abs+1)); % should I also adjust by nfreq?
                
                if strcmp(text_out,'') == 1
                    errorbar(times(tt),...
                    mean(n_curr(:,tt)),...
                    std(n_curr(:,tt)),...
                    'k','linewidth',4);
                    hold on
                else
                    errorbar(times(tt),...
                    mean(n_curr(:,tt)),...
                    std(n_curr(:,tt)),...
                    'g','linewidth',4);
                    hold on
                end
            end
            

            
            if t == 2 && f == 4
                 xlabel('Time relative to spike peak (s)')
            end 
           
            if t == 1 && strcmp(net_name,'coherence') == 1
                title(sprintf('%s',...
                    strrep(stats(n).time(t).freq(f).name,'_',' ')))
            elseif t == 1 && strcmp(net_name,'simple') == 1
                title('correlation')
            end
            
            xlim([-3 3]) 
            
        end
        
        
        
    end
end

for sp = 1:length(ha)
    axes(ha(sp))
    % formatting
    ylim([min_z-1 max_z+1]);
    if mod(sp,n_freq_abs+1) == 1
        ylabel(sprintf('%s s time window',stats(1).time(floor((sp/n_freq_abs+1))).name));
    else
        yticklabels([])
    end
    
    set(gca,'fontsize',15)
end

print(gcf,[out_folder,'nbs_change_',graph_method],'-depsc');

%% Say the patients with significant pre-spike rise
midpoint = nchunks/2;
for n = 1:length(network_names)
    fprintf('\n for %s:\n\n',network_names{n});
for t = 1:length(time_names)
    fprintf('\n for time window %s:\n\n',time_names{t});
    for i = 1:length(pt_names)
        fprintf('\n%s had significant pre-spike network change for:',pt_names{i});
        for f = 1:nfreq
            for tt = 1:midpoint - 1
            
                if i>size(stats(n).time(t).freq(f).p_all,1), continue; end
                % Get the p-value
                p = stats(n).time(t).freq(f).p_all(i,tt);
                
                if p < 0.05/(n_freq_abs+1)/length(pt_names)/(nchunks-1)
                    fprintf('\n%s time %d.\n',freq_names{f},tt);
                end

            end
        end
        fprintf('\n\n\n');
    end
end
end

end