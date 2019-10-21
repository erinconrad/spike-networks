function adj = get_simple_corr(values)
% This will just get a simple, bidirected pair-wise channel by channel
% correlation for each time signal. I am not doing partial correlation
% because I would like to compare time point to time point.


nchs = size(values,2);
adj = zeros(nchs,nchs);

for i = 1:nchs
    for j = 1:i-1
        c = abs(corr(values(:,i),values(:,j))); % absolute value of Pearson correlation
        adj(i,j) = c;
        adj(j,i) = c;
    end
end

end