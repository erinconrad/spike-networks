function out = get_clinical_info(orig)

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
pt_struct_folder = [main_folder,'data/spike_structures/'];

%% Load
pt = load([pt_struct_folder,'pt.mat']);
pt = pt.pt;

npts = length(orig.freq(1).pt);
for p = 1:npts
    name = orig.freq(1).pt(p).name;
    out(p).name = name;
    for i = 1:length(pt)
        if strcmp(pt(i).name,name)
            clinical = pt(i).clinical;
            out(p).clinical = clinical;
        end
    end
end

end