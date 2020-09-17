function plot_avg_spike(overwrite,not_a_spike)

%% Parameters
do_notch = 1; % notch filter?
do_car = 1; % common average reference?
pre_whiten = 0; % remove the AR(1) component for a pre-whitening step?

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
out_folder = [results_folder,'plots/avg_spikes/'];
sig_dev_folder = [results_folder,'signal_deviation/manual/'];


if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

% Load spike file for one patient to get the surround time
spike = load([eeg_folder,'HUP074_eeg.mat']);
spike = spike.spike;
surround_time = spike(1).surround_time;

pt = load(pt_file);
pt = pt.pt;

% get full directory listing
listing = dir(eeg_folder);
count = 0;
pt_names = {};

if not_a_spike
    not_text = '_not_spike';
else
    not_text = '';
end
    
for i = 1:length(listing)
    
    % Get name
    fname = listing(i).name;
    if contains(fname,'_eeg.mat') == 0, continue; end
    
    if not_a_spike
        if contains(fname,'not') == 0, continue; end
    else
        if contains(fname,'not') == 1, continue; end
    end
    pt_name = strsplit(fname,'_');
    pt_name = pt_name{1};
    
    if overwrite == 0
        if exist([out_folder,pt_name,not_text,'.eps'],'file') ~= 0
            fprintf('\nSkipping %s...\n',pt_name);
            continue
        else
            fprintf('\nDoing %s...\n',pt_name);
        end
    end
    
    pt_names = [pt_names;pt_name];
    
    % Find corresponding patient in pt struct
    pt_id = 0;
    for k = 1:length(pt)
        if strcmp(pt(k).name,pt_name) == 1
            pt_id = k;
            break
        end
    end
    if pt_id == 0
        fprintf('\nWarning, cannot find id for %s.\n',pt_name);
    end
    
    % Load the file
    spike = load([eeg_folder,fname]);
    spike = spike.spike;
    nspikes = length(spike);
    
    dev_all = zeros(size(spike(1).data,1),nspikes);
    dev_all_involved = zeros(size(spike(1).data,1),nspikes);
    dev_all_ch = zeros(size(spike(1).data,1),nspikes);
    surround = spike(1).surround_time;
    fs = spike(1).fs;
    
    ch_devs_all = zeros(nspikes,1);
    
    
    for s = 1:length(spike)
        
        % biggest deviation channel
        biggest_dev = spike(s).biggest_dev;
        all_involved = spike(s).involved;
        
        x = spike(s).data;
        
        % pre process
        x = pre_processing(x,do_car,pre_whiten,do_notch,fs);
        
        % all involved
        y = x(:,all_involved);
        
        % all channels
        z = x;
        
        % just most involved ch
        x = x(:,biggest_dev);
        
        
        
              
        % get baseline
        bl = median(x,1);
        bly = median(y,1);
        blz = median(z,1);
        
        % get dev
        dev = (abs(x-bl)).^2;
        dev_involved = (abs(y-bly)).^2;
        dev_involved = nanmean(dev_involved,2); % avg across involved chs
        dev_avg = (abs(z-blz)).^2;
        dev_avg = nanmean(dev_avg,2); % avg across all channels
                
        
        if size(dev_all,1) > size(dev,1)
            dev = [dev;...
                repmat(dev(end),size(dev_all,1)-size(dev,1),1)];
        elseif size(dev_all,1) < size(dev,1)
            dev(end-(size(dev,1)-size(dev_all,1))+1:end) = [];
        end
        
        if size(dev_all_involved,1) > size(dev_involved,1)
            dev_involved = [dev_involved;...
                repmat(dev_involved(end),size(dev_all_involved,1)-size(dev_involved,1),1)];
        elseif size(dev_all_involved,1) < size(dev_involved,1)
            dev_involved(end-(size(dev_involved,1)-size(dev_all_involved,1))+1:end) = [];
        end
        
        if size(dev_all_ch,1) > size(dev_avg,1)
            dev_avg = [dev_avg;...
                repmat(dev_avg(end),size(dev_all_ch,1)-size(dev_avg,1),1)];
        elseif size(dev_all_ch,1) < size(dev_avg,1)
            dev_avg(end-(size(dev_avg,1)-size(dev_all_ch,1))+1:end) = [];
        end
        
        dev_all(:,s) = dev;
        dev_all_involved(:,s) = dev_involved;
        dev_all_ch(:,s) = dev_avg;

        ch_devs_all(s) = biggest_dev;
        
    end
    
    final_dev = nanmean(dev_all,2);
    final_dev_involved = nanmean(dev_all_involved,2);
    final_dev_all = nanmean(dev_all_ch,2);
    
    % take mode channel with biggest dev
    ch_devs_all(ch_devs_all==0) = [];
    mode_ch = mode(ch_devs_all);

    soz_text = '';
    
    if pt_id ~= 0
        if isfield(pt(pt_id),'newSOZChs')
            if ~isempty(pt(pt_id).newSOZChs)
                soz = pt(pt_id).newSOZChs;
                
                if strcmp(spike(1).chLabels{mode_ch},...
                        pt(pt_id).new_elecs.names{mode_ch}) == 0
                    error('Ch names do not line up for %s\n',pt_name);
                end
                
                if ismember(mode_ch,soz) == 1
                    soz_text = 'is soz';
                else
                    soz_text = 'is not soz';
                end
            end
        end
    end
    
    
    
    figure
    set(gcf,'position',[440 58 932 740])
    subplot(2,1,1)
    plot(linspace(-surround,surround,length(final_dev)),final_dev)
    title(sprintf('%s avg spike power\nin most involved ch\nmode channel %s',pt_name,soz_text))
    set(gca,'fontsize',20)
    
    subplot(2,1,2)
    plot(linspace(-surround,surround,length(final_dev)),final_dev_all)
    title(sprintf('%s avg spike power\nin all chs',pt_name))
    set(gca,'fontsize',20)
    print([out_folder,pt_name,not_text],gcf,'-depsc');
    close(gcf)
        
end


end