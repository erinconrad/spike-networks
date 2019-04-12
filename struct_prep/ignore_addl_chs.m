function pt = ignore_addl_chs(pt)

%{
This is a one time use file to make sure I am ignoring the appropriate
electrodes and that coordinates are properly aligned
%}

[~,~,jsonfile,electrode_folder] = spike_network_files;

pt_info = loadjson(jsonfile);
ptnames = fieldnames(pt_info.PATIENTS);

%% Loop through patients in json file
for i = 1:length(ptnames)
    info = pt_info.PATIENTS.(ptnames{i});
    ignore_addl = info.IGNORE_ELECTRODES;
    
    % Find the appropriate pt and add list of ignored electrodes
    found = 0;
    for j = 1:length(pt)
        if strcmp(pt(j).name,ptnames{i}) == 1
            found = 1;
            pt(j).ignore.names = {};
            for k = 1:length(ignore_addl)
                pt(j).ignore.names = [pt(j).ignore.names;ignore_addl{k}];
            end
        end
    end
    
    if found == 0
        fprintf('Warning, did not find %s\n',ptnames{i});
    end
    
end


%% Loop through electrode data and reconstruct electrode info
for i = 1:length(pt)
    if isempty(pt(i).electrodeData) == 1, continue; end
    count = 0;
    pt(i).new_elecs.locs =[];
    pt(i).new_elecs.names = {};
    for j = 1:length(pt(i).electrodeData.electrodes)
        name = pt(i).electrodeData.electrodes(j).name;
        type = pt(i).electrodeData.electrodes(j).type;
        xyz = pt(i).electrodeData.electrodes(j).xyz;
        % don't include it in new elecs if it's in json ignore
        if ismember(name,pt(i).ignore.names) == 0
            count = count + 1;
            pt(i).new_elecs.electrodes(count).name = name;
            pt(i).new_elecs.electrodes(count).type = type;
            pt(i).new_elecs.electrodes(count).xyz = xyz;
            pt(i).new_elecs.locs = [pt(i).new_elecs.locs;xyz];
            pt(i).new_elecs.names = [pt(i).new_elecs.names;name];
        end
    end
    
    % Sanity checks
    %{
    if length(pt(i).ignore.names) ~= length(pt(i).electrodeData.electrodes) - ...
            length(pt(i).new_elecs.electrodes)
        error('what\n');
        
    end
    %}
    
    for j = 1:length(pt(i).ignore.names)
        if ismember(pt(i).ignore.names{j},pt(i).new_elecs.names) == 1
            error('what\n');
        end
    end
    
    for j = 1:length(pt(i).new_elecs.electrodes)
        name = pt(i).new_elecs.electrodes(j).name;
        for k = 1:length(pt(i).electrodeData.electrodes)
            if strcmp(name,pt(i).electrodeData.electrodes(k).name) == 1
                if isequal(pt(i).new_elecs.electrodes(j).xyz,...
                        pt(i).electrodeData.electrodes(k).xyz) == 0
                    error('what\n');
                end
            end
        end
    end
 
end

end