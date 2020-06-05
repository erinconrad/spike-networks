function sp = get_manual_times_from_excel

%{
This takes manual spike times from excel and fills a structure with spike
times
%}

%% Parameters
expected_num = 50; % expected number of spikes

%% Get locations
locations = spike_network_files;
script_folder = locations.script_folder;
main_folder = locations.main_folder;
addpath(genpath(script_folder));
sp_folder = [main_folder,'data/manual_spikes/'];
pt_folder = [main_folder,'data/spike_structures/'];

pt = load([pt_folder,'pt.mat']);
pt = pt.pt;


T = readtable([sp_folder,'manual spikes.xlsx']);

% Loop through variable names
for i = 1:length(T.Properties.VariableNames)
    name = T.Properties.VariableNames{i}; % pt name
    
    curr_pt = 0;
    % Get the pt number
    for whichPt = 1:length(pt)
        if strcmp(name,pt(whichPt).name) == 1
            curr_pt = whichPt;
            break
        end
    end
    
    if curr_pt == 0
        fprintf('Warning, could not find pt %s; skipping...\n',name); 
        continue
    end
    
    % Start filling in info for the spike structure
    sp(whichPt).name = name;
    
    % Get the spike time info
    sp_time_info = T.(name);
    
    % Figure out if there's text or if it's just an array
    type = class(sp_time_info);
    
    if strcmp(type,'double') == 1
        % It's just an array of spike times; add non-nans to array
        sp(whichPt).comment = [];
        sp(whichPt).spike = sp_time_info(~isnan(sp_time_info));
    elseif strcmp(type,'cell') == 1
        % There must be some text
        sp(whichPt).comment = {};
        sp(whichPt).spike = [];
        
        % Loop through each member of the cell array and add it to the
        % spike times if it is just a number
        for j = 1:length(sp_time_info)
            
            % Skip if empty
            if isempty(sp_time_info{j}) == 1, continue; end
            
            % See if it's a number
            if isempty(str2num(sp_time_info{j})) == 1
                % not a number; add to comments
                sp(whichPt).comment = [sp(whichPt).comment;sp_time_info{j}];
            else
                % number; add to spikes
                sp(whichPt).spike = [sp(whichPt).spike;str2num(sp_time_info{j})];
            end
        end
    else
        error('surprising input\n');
    end
    
    if length(sp(whichPt).spike) >= expected_num
        sp(whichPt).complete = 1;
    else
        sp(whichPt).complete = 0;
    end
end

if 0
    % show the spikes
    for i = 1:length(sp)
        if sp(i).complete == 1
            
            sp(i).name
            sp(i).comment
            sp(i).spike
            fprintf('\n\n');
        end
    end
end

end