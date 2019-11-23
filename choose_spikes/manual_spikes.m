function pt = manual_spikes

%{
Rules:
- no other spike 5 s before or 5 s after
%}


pt(3).name = 'HUP068';
pt(3).comment = 'looks like HFOs overriding spikes';
pt(3).spike(1).time = 741903.85;
pt(3).spike(2).time = 741979.49;
pt(3).spike(3).time = 742791.44;
pt(3).spike(4).time = 142404.32;
pt(3).spike(5).time = 143689.44;
pt(3).spike(6).time = 143962.63; 
pt(3).spike(7).time = 144061.47;
pt(3).spike(8).time = 211562.96;
pt(3).spike(9).time = 212233.02;
pt(3).spike(10).time = 256028.36;
%pt(3).spike(11).time = ;

end