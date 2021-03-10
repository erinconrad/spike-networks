function brain_and_elecs(pt,whichPt,ns)

name = pt(whichPt).name;
locs = pt(whichPt).new_elecs.locs;
offset = [-3 27 7.7653];

%% Get paths
gifti_script_path = '/Users/erinconrad/Desktop/research/spike_locations/scripts/tools/gifti-1.8';
addpath(genpath(gifti_script_path));

gifti_data_path = '/Users/erinconrad/Desktop/research/spike_locations/data/brains/';
giftiFolder = [gifti_data_path,name,'/'];
names = dir([giftiFolder,'*pial.gii']);
fname2 = names(1).name;
g = gifti([giftiFolder,fname2]);


transformed_elecs_path = '/Users/erinconrad/Desktop/research/spike_locations/data/transformedElectrodes/';
names = dir([transformed_elecs_path,name,'*']);
fname = names.name;
fileID = fopen([transformed_elecs_path,fname]);
out=textscan(fileID, '%s', 'whitespace',',');
out = out{1};
nchs = length(out)/5;

for i = 1:nchs
    electrodeData.electrodes(i).x = str2double(out{(i-1)*5+1});
    electrodeData.electrodes(i).y = str2double(out{(i-1)*5+2});
    electrodeData.electrodes(i).z = str2double(out{(i-1)*5+3});
    electrodeData.electrodes(i).xyz = [electrodeData.electrodes(i).x,...
       electrodeData.electrodes(i).y,electrodeData.electrodes(i).z];
    electrodeData.electrodes(i).name = out{(i-1)*5+4};
end


%% Align this with the original electrode data
% Number of electrodes different because my main electrode data ignores
% electrodes that have no ieeg data
newData(whichPt).locs = [];
for i = 1:length(pt(whichPt).new_elecs.electrodes)
    for j = 1:length(electrodeData.electrodes)
        if strcmp(pt(whichPt).new_elecs.electrodes(i).name,...
                electrodeData.electrodes(j).name) == 1
            newData(whichPt).electrodes(i).name = electrodeData.electrodes(j).name;
            newData(whichPt).electrodes(i).xyz = electrodeData.electrodes(j).xyz;
            newData(whichPt).locs = [newData(whichPt).locs;newData(whichPt).electrodes(i).xyz];
        end

    end

end

if length(pt(whichPt).new_elecs.electrodes) ~= length(newData(whichPt).electrodes)
    error('Number of electrodes not aligned\n');
end

newLocs = newData(whichPt).locs;
oldLocs = pt(whichPt).new_elecs.locs;
A = newLocs*pinv(oldLocs);

locs = A*locs-offset;
ns = ns';
p = plotGIFTI(g);
hold on
csz = 100;
scatter3(locs(:,1),locs(:,2),locs(:,3),csz,'k','linewidth',2);
scatter3(locs(:,1),locs(:,2),locs(:,3),csz,ns,'filled');
c = colorbar;
c.Label.String = 'Node strength';
c.FontSize = 20;
c.Position=[0.87 0.7033 0.0178 0.2467];
view(-120,-11);

%{
for k = 1:length(chs)
    scatter3(locs(k,1),locs(k,2),locs(k,3),350,this_col,'filled');
end
%}

end