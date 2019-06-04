function ch_nums = get_ch_nums_from_labels(pt,whichPt,ch_labels)

elecs = pt(whichPt).new_elecs.electrodes;
ch_nums = zeros(length(ch_labels),1);

for i = 1:length(ch_nums)
    for j = 1:length(elecs)
        if strcmp(ch_labels{i},elecs(j).name) == 1
            ch_nums(i) = j;
        end
    end
    
end

end