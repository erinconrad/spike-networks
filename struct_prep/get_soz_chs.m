function soz_chs = get_soz_chs(pt,name)

new_method = 1;

%{
if contains(name,'Study')
    new_method = 1;
else
    new_method = 2; % should give same result if 1 or 0; 2 means to use the new soz channels
end
%}


%% Find which pt
for i = 1:length(pt)
    if strcmp(pt(i).name,name)
        p = i;
    end
        
end

%fprintf('\nSOZ channels:\n');
if new_method == 1

    %% Get newsozchs (in old electrode space)
    soz_chs_old = pt(p).newSOZChs;

    %% Convert them into old electrode space labels
    soz_labels_old = cell(length(soz_chs_old),1);

    for i = 1:length(soz_chs_old)
        ch = soz_chs_old(i);
        % Get old electrode space label
        %fprintf('%s\n',pt(p).electrodeData.electrodes(ch).name);
        soz_labels_old{i} = pt(p).electrodeData.electrodes(ch).name;
    end

    %% Get the new indices of these labels (noting some will have none)
    soz_chs_new = nan(length(soz_chs_old),1);
    for i = 1:length(soz_labels_old)

        % Loop through new labels and find corresponding ch
        for j = 1:length(pt(p).new_elecs.electrodes)
            curr_name = pt(p).new_elecs.electrodes(j).name;

            if strcmp(soz_labels_old{i},curr_name)
                soz_chs_new(i) = j;
                break
            end
        end
    end

    num_nan = sum(isnan(soz_chs_new));
    soz_chs_new(isnan(soz_chs_new)) = [];
    if num_nan > 0
        fprintf('\nRemoved %d soz chs to fit new electrode structure for %s.\n',num_nan,name);
    end

    soz_chs = soz_chs_new;

elseif new_method == 0


    % Old approach
    %% Get sz labels
    sz_labels = {};

    for s = 1:length(pt(p).sz)
        for l = 1:length(pt(p).sz(s).electrodes)

            sz_labels = [sz_labels,pt(p).sz(s).electrodes{l}];
        end
    end

    sz_labels = unique(sz_labels);

    %% Get chs
    soz_chs = [];
    for s = 1:length(sz_labels)
        found_it = 0;
        for e = 1:length(pt(p).new_elecs.electrodes)
            if strcmp(pt(p).new_elecs.electrodes(e).name,sz_labels(s))
                soz_chs = [soz_chs,e];
                found_it = 1;
                break
            end
        end
       % if found_it == 0, error('what'); end
    end
    
elseif new_method == 2
    out_labels = turn_new_soz_into_struct(name);

    %% Get the new indices of these labels (noting some will have none)
    soz_chs = nan(length(out_labels),1);
    for i = 1:length(out_labels)

        % Loop through new labels and find corresponding ch
        for j = 1:length(pt(p).new_elecs.electrodes)
            curr_name = pt(p).new_elecs.electrodes(j).name;

            if strcmp(out_labels{i},curr_name)
                soz_chs(i) = j;
                break
            end
        end
    end


end

end