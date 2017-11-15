% Copied from "Matlab codes for peaks.docx" 01.11.17

% We're missing variables? 
% power1_100

% and functions?
% filter1_100
% findpeaks
% text
% movingstd2

thresh = 3; % cut off standard deviation
space = 10;
 
filter1_100=mean(power1_100)+thresh*std(power1_100);
subplot(3,5,1)
plot(freq1_100,power1_100)
[pks1,locs1]=findpeaks(power1_100,freq1_100,'MinPeakProminence',filter1_100,'Annotate','extents','MinPeakDistance',space);
text(locs1+.02,pks1,num2str((1:numel(pks1))'))

plot(freq200, power200)
[pks2,locs2]=findpeaks(power200,freq200,'MinPeakProminence',filter200,'Annotate','extents','MinPeakDistance',space);
text(locs2+.02,pks2,num2str((1:numel(pks2))'))



%% Code for Moving average and STDV:
M=movingstd2(power,20);
