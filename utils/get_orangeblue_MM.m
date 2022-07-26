function [ col ] = get_orangeblue_MM( n )
% [ col ] = get_orangeblue_MM( n )
% Create colormap for heatmaps, ranging from orange to blue
%
% Matteo Mischiati
% 2020

M=zeros(256,3);

%I will linearly interpolate between start and mid colors, and then between
%mid colors and end colors

start_red = 205; %230;
mid_red = 256;
end_red = 0;

M(1:128,1)=start_red + (1:128)*((mid_red-start_red)/128);
M(129:256,1)=mid_red + (1:128)*((end_red-mid_red)/128);

start_green = 128; %154;
mid_green = 256;
end_green = 128; %154;

M(1:128,2)=start_green + (1:128)*((mid_green-start_green)/128);
M(129:256,2)=mid_green + (1:128)*((end_green-mid_green)/128);

start_blue = 0;
mid_blue = 256;
end_blue = 205; %230;

M(1:128,3)=start_blue + (1:128)*((mid_blue-start_blue)/128);
M(129:256,3)=mid_blue + (1:128)*((end_blue-mid_blue)/128);

M = (1/256)*M;

for k = 1:n
    col(k,:) = M(round(1+(k-1)*255/(n-1)),:);
end

end