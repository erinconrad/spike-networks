function ers_change

%% Parameters

freq_bands = [0 4;... %delta
    4 8;...%theta
    8 12;...% alpha
    12 24;... %beta
    30 40;... % low gamma
    96 106;... % high gamma
    106 256;... %ultra-high
    0 256;... %broadband    
    ]; 
freq_names = {'delta','theta','alpha','beta','low_gamma',...
    'high_gamma','ultra_high','broadband'};
n_f = size(freq_bands,1);

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
eeg_folder = [main_folder,'results/eeg_data/'];
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
ers_folder = [results_folder,'ers/'];
out_folder = [results_folder,'plots/'];

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end


% Load spike file for one patient to get the surround time
spike = load([eeg_folder,'HUP074_eeg.mat']);
spike = spike.spike;
surround_time = spike(1).surround_time;


% Loop through time scales
n_freq_abs = 0;
max_F = 0;
network_count = 1;


% Loop through time scales
time_listing = dir(ers_folder);
time_count = 0;

for k = 1:length(time_listing)
    time_name= time_listing(k).name;
    time_window = str2num(time_name);

    % Skip if . or ..
    if strcmp(time_name,'.') == 1 || strcmp(time_name,'..') == 1
        continue
    end

    % Skip if not a directory
    if time_listing(k).isdir == 0, continue; end

    time_count = time_count + 1;
    stats(network_count).time(time_count).name = time_name;
    stats(network_count).time(time_count).time_window = time_window;
    time_folder = [ers_folder,time_name,'/'];

    pt_listing = dir([time_folder,'*.mat']);

    % load one to get nfreq
    ers = load([time_folder,pt_listing(1).name]);
    ers = ers.ers;
    nfreq = size(ers.freq_bands,1);
    if n_freq_abs < nfreq
        n_freq_abs = nfreq;
    end

    for f = 1:nfreq
        stats(network_count).time(time_count).freq(f).avg_ers_involved = ...
            nan(length(pt_listing),surround_time*2/time_window);

        stats(network_count).time(time_count).freq(f).t_all = ...
            nan(length(pt_listing),surround_time*2/time_window);
        
        stats(network_count).time(time_count).freq(f).p_all = ...
            nan(length(pt_listing),surround_time*2/time_window);
    end

    pt_names = {};
    
    % loop through pts
    for i = 1:length(pt_listing)

        pt_name = pt_listing(i).name;
        pt_name_pt = strsplit(pt_name,'_');
        pt_name_pt = pt_name_pt{1};
        pt_names = [pt_names;pt_name_pt];

        % load pt file
        ers = load([time_folder,pt_name]);
        ers = ers.ers;

        for f = 1:nfreq
            stats(network_count).time(time_count).freq(f).name = ers.freq_names{f};
            stats(network_count).time(time_count).freq(f).avg_ers_involved(i,:) = ers.powers_avg_involved(:,f);

            ers_involved = ers.powers_involved(:,:,f);
            
            % Do t-test
            for t = 2:size(ers_involved,2)
                [~,p,~,stats1] = ttest(ers_involved(:,1),ers_involved(:,t));
                stats(network_count).time(time_count).freq(f).t_all(i,t) = stats1.tstat;
                stats(network_count).time(time_count).freq(f).p_all(i,t) = p;
            end
            
            
            
        end


    end
    
    
    % combine t stats across patients to get group-level significance
    for f = 1:nfreq
        stats(network_count).time(time_count).freq(f).group_p = nan(surround_time*2/time_window,1);
        for t = 2:surround_time*2/time_window
            [~,p] = ttest(stats(network_count).time(time_count).freq(f).t_all(:,t));
            stats(network_count).time(time_count).freq(f).group_p(t) = p;
        end
        
        if 0
        fprintf('\nFor %s frequency:\n',freq_names{f});
       (array2table(stats(network_count).time(time_count).freq(f).p_all,...
            'RowNames',pt_names))
        end
    end
    %}
    
    %{
    for f = 1:nfreq
        stats(network_count).time(time_count).freq(f).group_p = nan(surround_time*2/time_window,1);
        for t = 2:surround_time*2/time_window
            
            [~,p] = ttest(stats(network_count).time(time_count).freq(f).avg_ers_involved(:,t),...
                stats(network_count).time(time_count).freq(f).avg_ers_involved(:,1));
         
            stats(network_count).time(time_count).freq(f).group_p(t) = p;
        end
    end
    %}

end
 
%% Plot
% Initialize figure
%{
nfreq columns and 2 (2 time scales) rows
%}
figure
set(gcf,'position',[100 100 1300 500])
[ha, pos] = tight_subplot(time_count, n_freq_abs, [0.01 0.01], [0.1 0.05], [0.05 0.01]);

for t = 1:time_count
    % change times for x axis
    nchunks = size(stats.time(t).freq(1).t_all,2);
    times = realign_times(nchunks,surround_time);
    
    nfreq = length(stats.time(t).freq);
    for f = 1:nfreq
        
        % this adds the number of frequencies + 1 if it's on the 2nd
        % time point (to move down a row), and it adds which frequency
        sp = (n_freq_abs)*(t-1) + f;
        axes(ha(sp));
        
        ers_curr = stats.time(t).freq(f).avg_ers_involved;
        % loop over patients and plot
        for i = 1:size(ers_curr,1)
            plot(times,ers_curr(i,:),'ko');
            hold on
        end
        
        % plot mean ers across patients
        for tt = 1:size(ers_curr,2)
            plot([times(tt)-0.25 times(tt)+0.25],...
                [mean(ers_curr(:,tt)) mean(ers_curr(:,tt))],...
                'k','linewidth',2);
        end
        
        % asterisks
        p = stats.time(t).freq(f).group_p;
        for tt = 2:length(p)
            text_out = get_asterisks(p(tt),(nchunks-1)*(n_freq_abs)); % should I also adjust by nfreq?
            text(times(tt),max(ers_curr(:,tt))+0.5,sprintf('%s',text_out),'fontsize',20,...
                    'horizontalalignment','center')
        end
        
        if t == 2 && f == 4
             xlabel('Time relative to spike peak (s)')
        end 
        
        title(sprintf('%s',...
                strrep(stats.time(t).freq(f).name,'_',' ')))

    end
    
end

%% Say the patients with significant pre-spike rise
midpoint = nchunks/2;
for t = 1:(time_count)
    fprintf('\n for time window %s:\n\n',time_listing(t).name);
    for i = 1:length(pt_names)
        fprintf('\n%s had significant pre-spike ERS for:',pt_names{i});
        for f = 1:nfreq
            for tt = 1:midpoint - 1
            
                % Get the p-value
                p = stats.time(t).freq(f).p_all(i,tt);
                
                if p < 0.05/nfreq/(nchunks-1)
                    fprintf('\n%s time %d.\n',freq_names{f},tt);
                end

            end
        end
        fprintf('\n\n\n');
    end
end


end