function view_utah

plot_time = [0 15]; % time to plot in seconds

%% Parameters
do_notch = 0; % notch filter?
do_car = 1; % common average reference?
pre_whiten = 0; % remove the AR(1) component for a pre-whitening step?

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
columbia_folder = locations.columbia_folder;
grid_folder = [columbia_folder,'Patient1_grid/'];
sleep_utah_folder = [grid_folder,'interictal_sleep_utah_0160_043647/'];

listing = dir(sleep_utah_folder);

for i = 1:length(listing)
    if contains(listing(i).name,'.mat') == 0, continue; end
    fname = listing(i).name;
    data = load([sleep_utah_folder,fname]);
    data = data.data;
    
    %% Grab basic info
    fs = data.MetaTags.SamplingFreq;
    elecs = data.ElectrodesInfo;
    
    %% EEG data
    values = double(data.Data);
    plot_idx = max([1,round(plot_time(1)*fs)]):...
        min([size(values,2),round(plot_time(2)*fs)]);
    nchs = size(values,1);
    
    %% Pre-processing
    temp_values = pre_processing(values',do_car,pre_whiten,do_notch,fs);
    values = temp_values'; % pre-processing assumes it's ntimepointsxnchs
    
    figure
    set(gcf,'position',[75 53 1311 752]);
    offset = 0;
    for ich = 1:nchs
        plot(linspace(plot_time(1),plot_time(2),length(plot_idx)),...
            values(ich,plot_idx)+offset);
        hold on
        if ich ~= nchs
            offset = offset - abs((prctile(values(ich+1,plot_idx),95) - prctile(values(ich,plot_idx),95)));
            %offset = offset - abs((max(values(ich+1,plot_idx)) - min(values(ich,plot_idx))));
        end
    
    end
    
    
    
end

end