function show_random_spikes(pt,cluster,whichPts)

%% Parameters
dur = 6;

%% Get locations
locations = spike_network_files;
script_folder = locations.script_folder;
main_folder = locations.main_folder;
pwname = locations.pwfile;
addpath(genpath(script_folder));

if isempty(locations.ieeg_folder) == 0
    addpath(genpath(locations.ieeg_folder));
end

%% Define which patients (all with >= 1 good cluster)
if isempty(whichPts)
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            if size(cluster(i).bad_cluster) < cluster(i).k
                whichPts = [whichPts,i];
            end
        end
    end
end

for whichPt = whichPts
    %% Get cluster info and remove bad clusters   
    all_spike_times = cluster(whichPt).all_times_all; % all spike times
    idx = cluster(whichPt).idx; % the cluster index for every spike
    all_spike_locs = cluster(whichPt).all_spikes; % spike locations
    bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad
    
    % Remove bad clusters
    bad_idx = (ismember(idx,bad_cluster));
    all_spike_times(bad_idx) = [];
    n_true_spikes = length(all_spike_times);
    all_spike_locs(bad_idx) = [];
    
    % Order spike times
    all_spike_times_old = all_spike_times;
    [all_spike_times,I] = sort(all_spike_times);
    all_spike_locs = all_spike_locs(I);
    
    while 1
        ix = randi(n_true_spikes);
        t = all_spike_times(ix);
        ch = all_spike_locs(ix);
        
        % Load the data
        fs = pt(whichPt).fs;
        indices = round((t-dur/2)*fs):round((t+dur/2)*fs);
        ieeg_name = pt(whichPt).ieeg_name;
        data = download_eeg(ieeg_name,indices,pwname,0,[]);
        
        chLabels = data.chLabels;
        values = data.values;
        nchs = size(chLabels,1);
        
       
        figure
        
        offset = 0;
        ch_offsets = zeros(nchs,1);
        ch_bl = zeros(nchs,1);

        ch_count = 0;
        for ich = 1:nchs
            ch_count = ch_count + 1;
            plot(linspace(0,dur,size(values,1)),values(:,ich)-offset,'k');
            ch_offsets(ch_count) = offset;
            ch_bl(ch_count) = -offset + median(values(:,ich));
            hold on
            text(dur+0.05,ch_bl(ch_count),sprintf('%s',chLabels{ich,1}))
            if ich<nchs
                if ~isnan(min(values(:,ich)) - max(values(:,ich+1)))
                    offset = offset - (min(values(:,ich)) - max(values(:,ich+1)));
                end
            end
        end

        pause
        close(gcf)
    end
    
end



end