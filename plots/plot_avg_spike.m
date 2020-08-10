function plot_avg_spike(overwrite)

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
for i = 1:length(listing)
    
    % Get name
    fname = listing(i).name;
    if contains(fname,'_eeg.mat') == 0, continue; end
    pt_name = strsplit(fname,'_');
    pt_name = pt_name{1};
    
    if overwrite == 0
        if exist([out_folder,pt_name,'.eps'],'file') ~= 0
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
    surround = spike(1).surround_time;
    fs = spike(1).fs;
    
    ch_devs_all = zeros(nspikes,1);
    
    for s = 1:length(spike)
        involved = spike(s).involved;
        
        if sum(involved) == 0, continue; end
        chs_involved = find(involved);
        
        x = spike(s).data;
        % pre process
        x = pre_processing(x,do_car,pre_whiten,do_notch,fs);
        
        
        x = x(:,involved);
      
        
        % get baseline
        bl = median(x,1);
        
        % get dev
        dev = abs(x-bl);
                
        % avg across spike chs
        avg_dev = mean(dev,2);
        
        if size(dev_all,1) > size(avg_dev,1)
            avg_dev = [avg_dev;...
                repmat(avg_dev(end),size(dev_all,1)-size(avg_dev,1),1)];
        elseif size(dev_all,1) < size(avg_dev,1)
            avg_dev(end-(size(avg_dev,1)-size(dev_all,1))+1:end) = [];
        end
        
        dev_all(:,s) = avg_dev;
        
        % Get max dev for each ch
        max_over_time = max(dev,[],1);
        
        % Find channel with biggest dev
        [~,biggest_dev_ch] = max(max_over_time);
        
        ch_devs_all(s) = chs_involved(biggest_dev_ch);
        
    end
    
    final_dev = mean(dev_all,2);
    
    % take mode channel with biggest dev
    ch_devs_all(ch_devs_all==0) = [];
    mode_ch = mode(ch_devs_all);
    
    mode_ch_dev = zeros(size(spike(1).data,1),nspikes);
    
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
                    soz_text = '(is soz)';
                else
                    soz_text = '(is not soz)';
                end
            end
        end
    end
    
    % get avg dev for that ch
    for s = 1:length(spike)
        x = spike(s).data;
        x = pre_processing(x,do_car,pre_whiten,do_notch,fs);
        x = x(:,mode_ch);
        
        bl = median(x,1);
        d = abs(x-bl);
        if size(mode_ch_dev,1) > size(d,1)
            d = [d;...
                repmat(d(end),size(mode_ch_dev,1)-size(d,1),1)];
        elseif size(d,1) > size(mode_ch_dev,1)
            d(end - (size(d,1)-size(mode_ch_dev,1))+1:end) = [];
        end
        mode_ch_dev(:,s) = d;
    end
    mode_ch_dev = mean(mode_ch_dev,2);
    
    figure
    set(gcf,'position',[440 58 932 740])
    subplot(2,1,1)
    plot(linspace(-surround,surround,length(final_dev)),final_dev)
    title(sprintf('%s',pt_name))
    set(gca,'fontsize',20)
    subplot(2,1,2)
    plot(linspace(-surround,surround,length(final_dev)),mode_ch_dev)
    title(sprintf('Biggest dev ch %s: %d',soz_text,mode_ch));
    set(gca,'fontsize',20)
    print([out_folder,pt_name],gcf,'-depsc');
    close(gcf)
        
end


end