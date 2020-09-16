function signal_dev_fig(windows,not_a_spike)

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

if not_a_spike
    not_a_spike_text = '_not_spike';
else
    not_a_spike_text = '';
end

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

% Load spike file for one patient to get the surround time
spike = load([eeg_folder,'HUP074_eeg.mat']);
spike = spike.spike;
surround_time = spike(1).surround_time;

% get full directory listing
listing = dir(sig_dev_folder);
count = 0;
for i = 1:length(listing)
    % look for only those that are directories, ignoring '.' and '..'
    if listing(i).isdir == 0
        continue;
    end
    
    if strcmp(listing(i).name,'.') == 1 || strcmp(listing(i).name,'..') == 1
        continue
    end
    
    
    
    % assume all other things are directories with different time windows
    
    time_text = listing(i).name;
    time_window = str2num(time_text);
    
    % Skip if not the time window we want
    if ismember(time_window,windows) == 0, continue; end
    
    time_window_folder = [sig_dev_folder,time_text,'/'];
    
    % load the file
    sub_listing = dir(time_window_folder);
    for k = 1:length(sub_listing)
        
        if contains(sub_listing(k).name,'.mat') == 0, continue; end
        
        if not_a_spike == 1
            if contains(sub_listing(k).name,'not_spike') == 0, continue; end
        else
            if contains(sub_listing(k).name,'not_spike') == 1, continue; end
        end
        count = count+1;
        temp_sig_dev = load([time_window_folder,sub_listing(k).name]);
        sig_dev(count).name = time_text;
        sig_dev(count).time_window = time_window;
        sig_dev(count).sig_dev = temp_sig_dev.sig_dev;
    end
end
n_windows = count;

%% Initialize arrays to get average z-scores
z_score_all = cell(length(sig_dev),1);
t_score_all = cell(length(sig_dev),1);
for k = 1:length(sig_dev)
    z_score_all{k} = zeros(length(sig_dev(k).sig_dev),...
    length(sig_dev(k).sig_dev(1).z_score_dev));

    t_score_all{k} = nan(length(sig_dev(k).sig_dev),...
    length(sig_dev(k).sig_dev(1).z_score_dev));
end
    

%% Initialize figure
figure
set(gcf,'Position',[300 500 1100 300]);
[ha, pos] = tight_subplot(1, n_windows, [0.05 0.05], [0.2 0.12], [0.05 0.05]);
for k = 1:n_windows
    axes(ha(k));
    if isfield(sig_dev(k).sig_dev(1),'time_window') == 1
        time_window = sig_dev(k).sig_dev(1).time_window;
        if length(time_window) >1
            nchunks = time_window;
            time_window = diff(time_window);
            time_window = time_window(1);
        else
            nchunks = size(sig_dev(k).sig_dev(1).index_windows,1);
        end
    else
        time_window = sig_dev(k).time_window;
        nchunks = size(sig_dev(k).sig_dev(1).index_windows,1);
    end
    time_text = sig_dev(k).name;
    curr_sig_dev = sig_dev(k).sig_dev;
    
   
    % change times for x axis
    times = realign_times(nchunks,surround_time);
    
    % Loop through patients and plot z scores at each time for each patient
    for i = 1:length(curr_sig_dev)
        z_score_temp = curr_sig_dev(i).z_score_dev;
        plot(times,z_score_temp,'ko');
        hold on
        
        % Add it to array
        z_score_all{k}(i,:) = z_score_temp;
        
        % Get t score and add this
        temp_temp_t = nan(length(curr_sig_dev(i).stats),1);
        for j = 2:length(curr_sig_dev(i).stats)
            temp_temp_t(j) = curr_sig_dev(i).stats(j).tstat;
        end
        t_score_all{k}(i,:) = temp_temp_t;
    end
    
    % plot the mean across patients
    curr_z_score = z_score_all{k};
    for j = 1:length(times)
        if j > 1
            [~,p] = ttest(t_score_all{k}(:,j)); % ttest looking at individual pt t stats
            %text_out = get_asterisks(p,size(z_score_all{k},2)-1);
            text_out = get_asterisks(p,1);
            if strcmp(text_out,'') == 1
                plot([times(j)-0.25*time_window times(j)+0.25*time_window],...
                [nanmean(curr_z_score(:,j)) nanmean(curr_z_score(:,j))],...
                'k','linewidth',4)
            else
                plot([times(j)-0.25*time_window times(j)+0.25*time_window],...
                [nanmean(curr_z_score(:,j)) nanmean(curr_z_score(:,j))],...
                'g','linewidth',4)
            end
        else
        
            plot([times(j)-0.25*time_window times(j)+0.25*time_window],...
                [nanmean(curr_z_score(:,j)) nanmean(curr_z_score(:,j))],...
                'k','linewidth',4)
        end
    end
    
    
   
    
    % formatting
    set(gca,'fontsize',20);
    xlabel('Time relative to spike peak (s)');
    if k == 1
        ylabel('Signal power (z-score)');
    end
    xlim([times(1)-0.25*time_window,times(end)+0.25*time_window]);
    title(sprintf('Time window %s s',time_text))
    ylim([min(min(z_score_all{k})) max(max(z_score_all{k}))+0.5])
    xl = get(gca,'xlim');
    xl(2) = 0.1;
    %xl(1) = -1;
    set(gca,'xlim',xl);
     %if k == 2, error('look\n'); end
end
    


print(gcf,[out_folder,'sig_dev',not_a_spike_text],'-depsc')

end