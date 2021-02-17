function is_depth = is_it_depth(pt,elec_num,name)

%% Find which pt
for i = 1:length(pt)
    if strcmp(pt(i).name,name)
        p = i;
    end
        
end

elec_name = pt(p).new_elecs.electrodes(elec_num).name;

if contains(elec_name,'D') &&...
        ~contains(elec_name,'GRID') &&...
        ~contains(elec_name,'DNET')
    is_depth = 1;
elseif contains(elec_name,'AMY') || contains(elec_name,'HIP')
    is_depth = 1;
elseif contains(elec_name,'LOF')
    is_depth = 1;
elseif contains(elec_name,'LA') && ...
        ~contains(elec_name,'LAT')
    is_depth = 1;
elseif contains(elec_name,'RA') && ...
        ~contains(elec_name,'RAT')
    is_depth = 1;
elseif contains(elec_name,'LH') || contains(elec_name,'RH')
    is_depth = 1;
else
    is_depth = 0;
end

end