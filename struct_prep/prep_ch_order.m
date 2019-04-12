function pt = prep_ch_order(pt)

% one time use function to get electrode order prior to getting ieeg data

[main_folder,pwfile,~,~,script_folder] = spike_network_files;
addpath(genpath(script_folder));

for whichPt = 1:length(pt)
    
    fprintf('Doing %s\n',pt(whichPt).name);
    
    % get ieeg name and sampling rate
    ieeg_name = pt(whichPt).ieeg_name;
    fs = pt(whichPt).fs;
    
    % Get fs and channel names
    data = download_eeg(ieeg_name,[],pwfile,1); 
    if data.fs ~= fs, error('What\n'); end

    % compare ch_labels to ones I want to ignore to get order to data I
    % want (so that it will match up with new_elecs)
    ch_labels = data.chLabels;
    ieeg_labels = ch_labels;
    for j = 1:length(ieeg_labels)
        ieeg_labels{j} = ieeg_ch_parser(ieeg_labels{j});
    end
    
    ch_order = [];
    label_order = {};
    
    for j = 1:length(pt(whichPt).new_elecs.electrodes)
        tname = pt(whichPt).new_elecs.electrodes(j).name;
        found = 0;
        for k = 1:length(ieeg_labels)
            if strcmp(ieeg_labels{k},tname) == 1
                found = 1;
                ch_order = [ch_order;k];
                label_order = [label_order;tname];
            end
        end
        if found == 0, error('what\n'); end
    end
    
    if isequal(label_order,pt(whichPt).new_elecs.names) == 0
        error('what\n');
    end
    
    pt(whichPt).new_elecs.ch_order = ch_order;
    
    
end