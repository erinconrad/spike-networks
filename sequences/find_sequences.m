function find_sequences


t2 = .05; % max time from preceding spike (15 ms in paper)
minSeqLength = 3; 

%% Parameters
surround_time = 2;
do_notch = 1; % notch filter?
do_car = 1; % common average reference?
pre_whiten = 0; % remove the AR(1) component for a pre-whitening step?

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
addpath(genpath(locations.BCT));
pt_file = [data_folder,'spike_structures/pt.mat'];
bct_folder = locations.BCT;
addpath(genpath(bct_folder));
% Folders
eeg_folder = [results_folder,'eeg_data/'];
seq_folder = [results_folder,'seq_data/'];
if exist(seq_folder,'dir') == 0
    mkdir(seq_folder)
end
listing = dir([eeg_folder,'*_eeg.mat']);

all_names = {};

for i = 1:length(listing)
    
    filename = listing(i).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    
    if contains(filename,'not') == 1, continue; end
    [a,b] = ismember(name,all_names);
    if a == 1
        pt_idx = b;
    else
        all_names = [all_names;name];
        pt_idx = length(all_names);
    end
    
    fprintf('\nDoing %s...\n',name);
    
    % load eeg data
    spike = load([eeg_folder,filename]);
    spike = spike.spike;
    values = spike(1).data;
    nchs = size(values,2);
    fs = spike(1).fs;
    t2_idx = fs*t2;
    seq = [];
    
    % loop through spikes
    for s = 1:length(spike)
        
        seq(s).seq = [];
        
                
        % get ordered chs
        order = spike(s).ordered_chs;
        
        times = order(:,2);
        chs = order(:,1);
        
        % Skip it if it's less than the minimum sequence length
        if size(order,1) < minSeqLength
            continue;
        end
        
        % Get the difference in time
        time_diff = diff(times);
        
        % Find those under the time limit
        under_limit = time_diff < t2;
        
        longest_span = 0;
        longest_span_idx = 0;
        % Take the longest span of ones
        
        % Loop through starting index
        for j = 1:length(under_limit)
            
            % find how long the span is starting with that index
            span = 0;
            for k = j:length(under_limit)
                if under_limit(k) == 1
                    span = span + 1;
                else
                    break % break at first non one
                end
            end
            
            if span > longest_span
                longest_span = span;
                longest_span_idx = j;
            end
        end
        
        if longest_span < minSeqLength
            continue;
        end
        
        curr_seq = order(longest_span_idx:longest_span_idx+longest_span-1,:);
        out_chs = curr_seq(:,1);
        out_times = curr_seq(:,2);
        
        % add deviation
        dev = spike(s).rel_dev(out_chs);
        seq(s).seq = [curr_seq,dev];
        
        
        % Add some superlatives
        [~,biggest_idx] = max(dev);
        biggest = out_chs(biggest_idx);
        seq(s).biggest_ch = biggest;
        
        % find all first chs
        first_time = out_times(1);
        first_idx = (out_times == first_time);
        % take the biggest as the tie-breaker
        size_ties = dev(first_idx);
        [~,tie_breaker] = max(size_ties);
        first = out_chs(tie_breaker);
        seq(s).first_ch = first;
        
        % Double check
        final_times = curr_seq(:,2);
        if any(diff(final_times>t2))
            error('what');
        end
        
        if 0
        data = spike(s).data;
        data = pre_processing(data,do_car,pre_whiten,do_notch,fs);
        figure
        set(gcf,'position',[100 100 1000 500])
        offset = 0;
        ch_offsets = zeros(size(curr_seq,1),1);
        ch_bl = zeros(size(curr_seq,1),1);
        for ich = 1:size(curr_seq,1)
            
            if out_chs(ich) == first
                plot(linspace(-3,3,size(data,1)),data(:,out_chs(ich))+offset,'g')
            elseif out_chs(ich) == biggest
                plot(linspace(-3,3,size(data,1)),data(:,out_chs(ich))+offset,'r')
            else
                plot(linspace(-3,3,size(data,1)),data(:,out_chs(ich))+offset,'k')
            end
            
            ch_offsets(ich) = offset;
            ch_bl(ich) = offset + median(data(:,out_chs(ich)));
            hold on
            text(surround_time+0.05,ch_bl(ich),sprintf('%s',spike(s).chLabels{out_chs(ich)}))        
            if ich<size(curr_seq,1)
                offset = offset + max(data(:,out_chs(ich))) - min(data(:,out_chs(ich+1)));
            end

        end
        xlim([-surround_time,surround_time]);
        title(sprintf('%s spike %d first ch %s',name,s,spike(s).chLabels{out_chs(1)}))
        pause
        close(gcf)
        end
        
    end
    
    % Save the file
    save([seq_folder,name,'_seq.mat'],'seq');
    
end

end