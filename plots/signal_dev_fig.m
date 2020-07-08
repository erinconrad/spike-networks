function signal_dev_fig

%% Parameters
alpha = 0.05;

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


if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end



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
    count = count+1;
    time_text = listing(i).name;
    time_window = str2num(time_text);
    time_window_folder = [sig_dev_folder,time_text,'/'];
    
    % load the file
    temp_sig_dev = load([time_window_folder,'sig_dev.mat']);
    sig_dev(count).name = time_text;
    sig_dev(count).time_window = time_window;
    sig_dev(count).sig_dev = temp_sig_dev.sig_dev;
    
end
n_windows = count;

%% Initialize arrays to get average z-scores
z_score_all = cell(length(sig_dev),1);
for k = 1:length(sig_dev)
    z_score_all{k} = zeros(length(sig_dev(k).sig_dev),...
    length(sig_dev(k).sig_dev(1).z_score_dev));
end
    

%% Initialize figure
figure
set(gcf,'Position',[300 500 1100 300]);
[ha, pos] = tight_subplot(1, n_windows, [0.05 0.05], [0.2 0.1], [0.1 0.05]);
for k = 1:n_windows
    axes(ha(k));
    time_window = sig_dev(k).time_window;
    time_text = sig_dev(k).name;
    curr_sig_dev = sig_dev(k).sig_dev;
    
    % Loop through patients and plot z scores at each time for each patient
    for i = 1:length(curr_sig_dev)
        z_score_temp = curr_sig_dev(i).z_score_dev;
        plot(z_score_temp,'ko');
        hold on
        
        % Add it to array
        z_score_all{k}(i,:) = z_score_temp;
    end
    
    % plot the mean across patients
    for j = 1:size(z_score_all{k},2)
        curr_z_score = z_score_all{k};
        plot([j-0.25 j+0.25],...
            [nanmean(curr_z_score(:,j)) nanmean(curr_z_score(:,j))],...
            'k','linewidth',2)
    end
    
    yl = get(gca,'ylim');
    
    % Do a t-test comparing the z-scores between the first and subsequent
    % time points
    for j = 2:size(z_score_all{k},2)
        [~,p] = ttest(z_score_all{k}(:,1),z_score_all{k}(:,j));
        text_out = get_asterisks(p,size(z_score_all{k},2));
        text(j,yl(2)-0.2,text_out,'fontsize',20,...
            'HorizontalAlignment','Center')
    end
    
end
    




end