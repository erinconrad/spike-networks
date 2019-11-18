function plot_avg_spike(whichPts)

%% Parameters
time_window = 0.5; %in seconds
n_chunks = 22; % number of time windows


%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
spike_times_file = [data_folder,'spike_times/times.mat'];
pt_file = [data_folder,'spike_structures/pt.mat'];

times = load(spike_times_file); % will result in a structure called "out"
times = times.out;
pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

if isempty(whichPts) == 1
    whichPts = 1:length(times);
end

 % output folder
out_folder = [results_folder,'plots/avg_dev/'];

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

% Loop through patients
for whichPt = whichPts
    
    %% Prep patient
    
    % Skip it if name is empty
    if isempty(times(whichPt).name) == 1, continue; end
    name = times(whichPt).name;
    fprintf('\nDoing %s\n',name);
    
    
    % Get basic data about patient
    pt_folder = [results_folder,name,'/'];
    data = load([pt_folder,'basic_info.mat']); % returns a structure called data
    data = data.data;
    
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
    
    
    
    % Loop through spike files
    s_count = 0;
    
    listing = dir([pt_folder,'spikes_*.mat']);
    for f = 1:length(listing)
        
        
        spike = load([pt_folder,listing(f).name]);
        spike = spike.spike;
 
        % Loop through spikes
        for s = 1:length(spike)
            %fprintf('Doing file %d of %d...\n',f,10);
            %fprintf('Doing spike %d of %d...\n',s,length(spike));
            %tic
            if isempty(spike(s).time) == 1, continue; end
            
            s_count = s_count + 1;
            
            % Grab the appropriate channels
            values = spike(s).values(:,~ignore);
            is_sp_ch = strcmp(ch_labels(~ignore),spike(s).label);
            
            %% Figure out times for which I will be calculating adjacencies
            % The peak should be the very center of each file
            peak = round(size(values,1)/2);
            
            % Get index windows
            index_windows = zeros(n_chunks,2);
            tick_window = time_window*data.fs;
            
            for i = 1:n_chunks
                index_windows(i,1) = peak - tick_window*n_chunks/2 + tick_window*(i-1);
                index_windows(i,2) = peak - tick_window*n_chunks/2 + tick_window*(i);
            end
            
            % For first spike, initialize array
            if f==1 && s == 1
                all_dev = zeros(1000,round(index_windows(end,2) - index_windows(1,1) +1));
                times_plot = [(index_windows(1,1)-peak)/data.fs,(index_windows(end,2)-peak)/data.fs];
            end
            
            % Get values from spike channel
            curr_values = values(index_windows(1,1):index_windows(end,2),is_sp_ch);
            
            % Get deviation from median for each point
            curr_dev = sqrt((curr_values - median(curr_values)).^2);
            
            % Fill array with deviation
            all_dev(s_count,:) = curr_dev;
            
        end
        
    end
    
    %% Plot average deviation across all spikes
    mean_dev = nanmean(all_dev,1);
    figure
    set(gcf,'position',[1 352 1375 446])
    plot(linspace(times_plot(1),times_plot(2),length(mean_dev)),mean_dev,'k','linewidth',2);
    hold on
    yl = get(gca,'ylim');
    for i = 1:n_chunks
        plot(([index_windows(i,1) index_windows(i,1)]-peak)/data.fs,[yl(1) yl(2)],'k--');
    end
    plot(([index_windows(end,2) index_windows(end,2)]-peak)/data.fs,[yl(1) yl(2)],'k--');
    yticklabels([])
    xlabel('Time relative to spike peak (s)')
    title(sprintf('%s average signal deviation from baseline',name))
    set(gca,'fontsize',20)
    print(gcf,[out_folder,name],'-depsc')
    close(gcf)
    
end