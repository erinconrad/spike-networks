function shift = find_subarray(A,b)

shift = nan;
for i = 1:length(A) - length(b) + 1
    if isequal(A(i:i+length(b)-1),b)
        shift = i - 1;
        break
    end
end

end