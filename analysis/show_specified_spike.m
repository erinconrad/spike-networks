function show_specified_spike(p,old_spike_indices,all_names,bad)

%{
Currently only works for spikes (not "not spikes", because I don't have the
not spike data on my laptop)
%}

figure
set(gcf,'position',[400  800  500 1000])

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
eeg_folder = [results_folder,'eeg_data/'];

for i = 1:length(old_spike_indices)
    old_s = old_spike_indices(i);

    %% Convert cleaned spike number to original spike number
    s = clean_to_orig_index(p,old_s,1,bad);

    %% get eeg data
    name = all_names{p};
    spike = load([eeg_folder,name,'_eeg.mat']);
    spike = spike.spike;
    values = spike(s).data;
    nchs = size(values,2);
    chLabels = spike(s).chLabels;

    %% Plot the eeg
    
    offset = 0;
    ch_offsets = zeros(nchs,1);
    ch_bl = zeros(nchs,1);

    ch_count = 0;
    for ich = 1:nchs
        ch_count = ch_count + 1;
        plot(linspace(0,6,size(values,1)),values(:,ich)-offset,'k');
        ch_offsets(ch_count) = offset;
        ch_bl(ch_count) = -offset + median(values(:,ich));
        hold on
        text(6+0.05,ch_bl(ch_count),sprintf('%s',chLabels{ich}))
        if ich<nchs
            if ~isnan(min(values(:,ich)) - max(values(:,ich+1)))
                offset = offset - (min(values(:,ich)) - max(values(:,ich+1)));
            end
        end
    end

    pause
    hold off
end

end
