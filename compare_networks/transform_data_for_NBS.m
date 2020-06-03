function [mat,design] = transform_data_for_NBS(meta,which_freq,which_time)
    
%{
This takes adjacency matrix data in my original structure and outputs
matrices in a format that can be run by NBS Connectome. 

meta is the structure with the adjacency matrices
which_freq is which frequency band
which_time is the comparison time relative to the first time

mat is an nch x nch x nspikes*2 matrix, where nch is the number of
channels, and each of the nspikes*2 matrices is the adjacency matrix for
that particular spike. There are nspike*2 matrices because there will be
both the first time and the comparison time.

design is an nspikes*2 x 2 matrix of 1s and 0s, where 1s in the first
column indicate that it is the first time and 1s in the second column
indicate that it is the second time.
%} 


nspikes = length(meta.spike);
nchs = size(meta.spike(1).adj(which_freq).adj,2);

mat = zeros(nchs,nchs,nspikes*2);
design = zeros(nspikes*2,2);

for t = 1:2
    for s = 1:nspikes
        all_times_adj = meta.spike(s).adj(which_freq).adj;
        if t == 1
            adj = squeeze(all_times_adj(1,:,:));
            design((t-1)*nspikes+s,1) = 1;
            design((t-1)*nspikes+s,2) = 0;
        else
            adj = squeeze(all_times_adj(which_time,:,:));
            design((t-1)*nspikes+s,1) = 0;
            design((t-1)*nspikes+s,2) = 1;
        end
        mat(:,:,(t-1)*nspikes+s) = adj;
        
    end
end

end