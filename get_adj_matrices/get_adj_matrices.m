function adj = get_adj_matrices(values,fs,freq_bands)

nchs = size(values,2);

% Initialize adjacency matrix
for i = 1:size(freq_bands,1)
    adj(i).adj = zeros(nchs,nchs);
end


% Do it for each frequency 
for i = 1:nchs
    for j = 1:i-1
        x = values(:,i);
        y = values(:,j);

       
        % This is the same as mscohere, except it removes an evalin call
        % that was taking up the majority of the computational time. I
        % expect that I cannot call this with additional arguments because
        % of this change.
        cxy = mscohere_erin(x,y);
        f = 0:2:fs/2;
        
      %  [cxy_old,f] = mscohere(x,y,[],[],[],fs);
        
        for ff = 1:size(freq_bands,1)
            adj(ff).adj(i,j) = mean(cxy(f>=freq_bands(ff,1)&...
                f<=freq_bands(ff,2)));
            adj(ff).adj(j,i) = mean(cxy(f>=freq_bands(ff,1)&...
                f<=freq_bands(ff,2)));
            
        end
        %}
    end
end

if 0
    figure
    subplot(1,2,1)
    imagesc(adj(1).adj)
    title('High gamma')
    colorbar
    
    subplot(1,2,2)
    imagesc(adj(2).adj)
    title('Beta')
    colorbar
    pause
    close(gcf)
end

end