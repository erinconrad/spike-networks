function manual_spikes

%{
Rules:
- no other spike 5 s before or 5 s after
%}

locations = spike_network_files;
main_folder = locations.main_folder;

out_folder = [main_folder,'data/manual_spikes/'];

%% HUP068

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

%% HUP070 
% I can't find any isolated spikes

%% HUO073

%% HUP074
pt(6).name = 'HUP074';
s = [449468.43
449513.24
449891.1
449939.85
213267.6
231614.88
464147.59
464516.62
464867.73
561522.75
562433.78
562703.17
563130.2];
for i = 1:length(s)
    pt(6).spike(i).time = s(i);
end

%% HUP075
whichPt = 7;
pt(whichPt).name = 'HUP075';
s = [761002.47
761024.13
761397.78
761432.89
762137.28
762670.6
595214.33
595308.68
595346.61
595408.99
596119.91
832943.29
832971.72];
for i = 1:length(s)
    pt(whichPt).spike(i).time = s(i);
end

%% HUP078
whichPt = 8;
pt(whichPt).name = 'HUP078';
s = [367617.45
367850.31
367916.97
367959.24
458492.75
458543.06
458550.01
458585.74
459160.22
459193.21
459224.57
459326.12
459737.78
459745.68];
for i = 1:length(s)
    pt(whichPt).spike(i).time = s(i);
end

%% HUP080
whichPt = 9;
pt(whichPt).name = 'HUP080';
s = [1109244.66
1109468.25
1109566.55
1109876.21
1109935.34
1109963.48
1109993.66
814833.15
814876.66
868740.74
1528770.81
1528872.09
1529003.43
1529015.53];
for i = 1:length(s)
    pt(whichPt).spike(i).time = s(i);
end

%% HUP087
whichPt = 13;
pt(whichPt).name = 'HUP087';
s = [125295.11
125377.37
125818.57
132250.15
132584.62
132670.83
132898.57
134617.74
151682.91
153560.06
154125.69
154289.27
203350.26];
for i = 1:length(s)
    pt(whichPt).spike(i).time = s(i);
end

%% HUP0105
whichPt = 16;
pt(whichPt).name = 'HUP105';
s = [943927.1
944152.59
944230.65
944522.96
433081.92
813247.6
170836.36
170901.22
171041.15
171303.22
171846.93
171878.39
235094.22
235247.39];
for i = 1:length(s)
    pt(whichPt).spike(i).time = s(i);
end

%% HUP0106
whichPt = 17;
pt(whichPt).name = 'HUP106';
s = [658756.02
659465.67
659691.09
659905.42
659932.1
557901.75
558300.73
558500.16
876417.38
876417.38
876585.55
394915.37
666561.75
667046.11];
for i = 1:length(s)
    pt(whichPt).spike(i).time = s(i);
end


%% Output
sp = pt;
save([out_folder,'sp.mat'],'sp');



end