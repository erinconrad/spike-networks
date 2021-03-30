function out = get_data_from_spec_ch(data,which_chs)

if length(which_chs) ~= size(data,1)
    error('Data must be of size n_spikes x n_times x n_chs');
end

if max(which_chs) > size(data,3)
    error('Data must be of size n_spikes x n_times x n_chs');
end

out = nan(size(data,1),size(data,2));

for s = 1:size(data,1)
    % get the data on the channel specified by which_ch
    out(s,:) = squeeze(data(s,:,which_chs(s)));
end

end