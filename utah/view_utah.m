function view_utah

%{
spikes on AST 4-5 and MST 4-5 (subtemporal strips)

utah array under temporal grid, between MTG15 and MTG16

I don't see spikes by MTG 15 or 16, and so I don't think we would expect to
see epileptiform discharges on the Utah array!
%}

plot_time = 15; % time to plot in seconds
goal_fs = 500;

%% Parameters
do_notch = 1; % notch filter?
do_car = 1; % common average reference?
pre_whiten = 0; % remove the AR(1) component for a pre-whitening step?

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
script_folder = locations.script_folder;
addpath(genpath(script_folder));
columbia_folder = locations.columbia_folder;
grid_folder = [columbia_folder,'Patient1_grid/'];
sleep_utah_folder = [grid_folder,'interictal_sleep_utah_0160_043647/'];
sleep_macro_folder = [grid_folder,'interictal_sleep_macro/'];

listing = dir(sleep_utah_folder);
%listing = dir(sleep_macro_folder);

figure
set(gcf,'position',[75 53 1311 752]);

for i = 1:length(listing)
    if contains(listing(i).name,'.mat') == 0, continue; end
    fname = listing(i).name;
    data = load([listing(i).folder,'/',fname]);
    struct_name = fieldnames(data);
    struct_name = struct_name{1};
    data = data.(struct_name);
    
    %% Grab basic info
    old_fs = data.MetaTags.SamplingFreq;
    elecs = data.ElectrodesInfo;
    
    %% EEG data
    values = double(data.Data);
    values = values';
    total_time = size(values,1)/old_fs;
    
    %% Downsample to 500 Hz
    n_downsample = round(old_fs/goal_fs);
    values = downsample(values,n_downsample);
    fs = size(values,1)/total_time;
    
    
    nchs = size(values,2);
    
    %% Pre-processing
    values = pre_processing(values,do_car,pre_whiten,do_notch,fs);

    
    
    start_idx = 1;
    
    while 1
        plot_idx = start_idx:...
            min([size(values,1),round(plot_time*fs)+start_idx]);
        

        offset = 0;
        for ich = 1:nchs
            plot(linspace(0,plot_time,length(plot_idx)),...
                values(plot_idx,ich)+offset);
            hold on
            text(plot_time + 0.2,median(values(plot_idx,ich)+offset),...
                sprintf('%s',elecs(ich).Label))
            if ich ~= nchs
                offset = offset - abs((prctile(values(plot_idx,ich+1),98) - prctile(values(plot_idx,ich),2)));
                %offset = offset - abs((max(values(ich+1,plot_idx)) - min(values(ich,plot_idx))));
            end

        end
        
        pause
        
        if size(values,1)>round(plot_time*fs)+start_idx
            start_idx = round(plot_time*fs)+start_idx;
            hold off
        else
            hold off
            break
        end
    
    end
    
    
    
end

end
