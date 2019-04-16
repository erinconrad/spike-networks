function spike_validator(whichPts)

%{
This function takes a bunch of downloaded eeg data and spike times,
randomly chooses 50 spikes, and plots the eeg data for the spike channel
and 3 randomly chosen other channels to get a sense of what the data looks
like

%}

%% Parameters
redo = 1; % 1 if we want to re-do the ones already done
n_spikes = 50; % number of spikes to plot
n_channels = 3; % number of channels besides spike channel to plot
n_chunks = 5; % How many plots per patient

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
spike_times_file = [data_folder,'spike_times/times.mat'];
pt_file = [data_folder,'spike_structures/pt.mat'];

times = load(spike_times_file); % will result in a structure called "times"
times = times.times;
pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

if isempty(whichPts) == 1
    whichPts = 1:length(times);
end

%% Loop through patients
for whichPt = whichPts
    
    % Skip it if name is empty
    if isempty(times(whichPt).name) == 1, continue; end
    name = times(whichPt).name;
    fprintf('Doing %s\n',name);
    
    % Make validation folder within pt folder
    pt_folder = [results_folder,name,'/'];
    val_folder = [pt_folder,'validation/'];
    if exist(val_folder,'dir') == 0, mkdir(val_folder); end
    
    if redo == 0
        listing = dir([val_folder,'*.eps']);
        if length(listing) == 6
            fprintf('Already did %s, skipping...\n',name);
            continue;
        end
    end
    
    % Select 50 random spikes
    n_spikes_total = sum(~isnan(times(whichPt).spike_times));
    which_spikes = sort(randi(n_spikes_total,n_spikes,1));
    
    % initialize out structure
    out = [];
    count = 0;
    
    % Get basic data about patient
    data = load([pt_folder,'basic_info.mat']); % returns a structure called data
    data = data.data;
    all_chs = 1:length(data.chLabels);
    
    % convert ch labels to nice form and decide which to ignore
    ch_labels = data.chLabels(:,1);
    ignore = zeros(length(ch_labels),1);
    
    for i = 1:length(ch_labels)
        ch_labels{i} = ieeg_ch_parser(ch_labels{i});
        for j = 1:length(pt(whichPt).ignore.names)
            if strcmp(pt(whichPt).ignore.names(j),ch_labels{i}) == 1
                ignore(i) = 1;
            end
        end
    end
    no_ignore = find(ignore == 0);
    
    % Select what files they are coming from
    which_files = floor((which_spikes-1)/100)+1;
    unique_files = unique(which_files);
    
    % Loop through the files (doing 1 file at a time since they're big so
    % they take a little while to load. They will be cleared each time I 
    % load a new one)
    for f = 1:length(unique_files)
        file_temp = sprintf('spikes_%d.mat',unique_files(f));
        spike = load([pt_folder,file_temp]); % this will create a structure called "spike"
        spike = spike.spike;
        
        % Identify the spike(s) associated with that file
        curr_spikes = (which_files == unique_files(f));
        curr_spikes = which_spikes(curr_spikes);
        
        % Loop through spikes and select appropriate channels
        for s = 1:length(curr_spikes)
            sp_idx = rem(curr_spikes(s),100);
            if sp_idx == 0, sp_idx = 100; end
            
            all_chs = 1:size(spike(sp_idx).values,2);
            
            % select 3 random channels other than the channel the spike is
            % on
            sp_ch = find(strcmp(ch_labels,spike(sp_idx).label));
            no_ignore(no_ignore==sp_ch) = [];
            rand_chs = randi(length(no_ignore),n_channels,1);
            rand_chs = no_ignore(rand_chs);
            all_chs = [sp_ch;rand_chs];
            
            % grab the data
            count = count + 1;
            out(count).data = spike(sp_idx).values(:,all_chs);
            out(count).labels = ch_labels(all_chs);
            if strcmp(out(count).labels{1},spike(sp_idx).label) ~= 1, error('what\n'); end
            out(count).time = spike(sp_idx).time;
            
        end
        
    end
    
    % Now that I have all 50 spikes, plot 5 at a time
    count = 0;
    for i = 1:n_chunks
        to_do = count+1:count+n_spikes/n_chunks;
        figure
        set(gcf,'position',[72 21 1300 780]);
        [ha,pos] = tight_subplot(n_spikes/n_chunks,1,[0.01 0],[0.01 0.01],[0.01 0.01]);
        for j = 1:length(to_do)
            idx = to_do(j);
            axes(ha(j));
            
            plot(out(idx).data(:,1),'k','linewidth',2);
            hold on
            yl = get(gca,'ylim');
            xl = get(gca,'xlim');
            plot([size(out(idx).data,1)/2 size(out(idx).data,1)/2],[yl(1) yl(2)],'k--');
            legend(out(idx).labels{1});
            
            text((xl(1)+xl(2))/2,(yl(1)+0.9*(yl(2)-yl(1))),sprintf('%1.1f s, window = 10 s',...
                out(idx).time-5),'fontsize',15);
            xticklabels([])
            yticklabels([])
        end
        
        fig_file_name = sprintf('Fig%d',i);
        print([val_folder,fig_file_name],'-depsc');
        close(gcf);
        count = count + 10;
    end
    
    % Make an additional plot showing the average deviation over time
    dev = zeros(n_spikes,size(out(1).data,1));
    for i = 1:length(out)
        t = sqrt((out(i).data(:,1)-median(out(i).data(:,1))).^2)';   
        if size(t,2) > size(dev,2) %sometimes it's 5001 points instead of 5000
            t = t(1:size(dev,2));
        end
        dev(i,:) = t;
    end
    avg_dev = nanmean(dev,1);
    figure
    set(gcf,'position',[72 21 1300 100]);
    plot(avg_dev,'k','linewidth',2)
    print([val_folder,'avg_dev'],'-depsc');
    pause(0.5)
    close(gcf);
    
end

end