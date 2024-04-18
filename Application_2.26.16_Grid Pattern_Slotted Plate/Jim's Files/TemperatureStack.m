function [Tstack] = TemperatureStack()

[infname, directory] = uigetfile('*.lst', 'Select List File');
filename = fullfile(directory,infname);

[ ImageNames ] = FileList( filename );
num = length(ImageNames);

in{1} = fullfile(directory,'background_Avg.tif');
in{2} = fullfile(directory,'wind-off_room temp_Avg.tif');
in{3} = fullfile(directory,ImageNames{1});
T0 = Point(in);
figure(1)
imagesc(T0,[40 80]);
colormap gray;
impixelinfo;


for i=2:num,
    in{3} = fullfile(directory,ImageNames{i});
    T = Point(in);
    dT = T-T0;
    Tstack(:,:,i-1) = dT;
    figure(2)
    imagesc(dT,[-10 1]);
    colormap gray;
    pause(1)
end


% File list function
function [ out ] = FileList( filename )
num = 1;
fid=fopen(filename);
while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end
    out{num} = tline;
    num = num+1;
end
fclose(fid);