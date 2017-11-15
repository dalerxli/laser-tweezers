% Copied from "Matlab codes for peaks.docx" 01.11.17

thresh = 3; % cut off standard deviation
space = 10;
 
%% introduction of power100 -- is this just normal PSD output from fftz.py?
filter100=mean(power100)+thresh*std(power100);
subplot(3,5,1)
plot(freq100,power100)
%% consider oct2py.octave, https://blog.ytotech.com/2015/11/01/findpeaks-in-python/ 
%% also: PeakUtils (peakutils.indexes), 
[pks1,locs1]=findpeaks(power100,freq100,'MinPeakProminence',filter100,'Annotate','extents','MinPeakDistance',space);
text(locs1+.02,pks1,num2str((1:numel(pks1))'))

plot(freq200, power200)
[pks2,locs2]=findpeaks(power200,freq200,'MinPeakProminence',filter200,'Annotate','extents','MinPeakDistance',space);
text(locs2+.02,pks2,num2str((1:numel(pks2))'))


%% Code for Moving average and STDV:
M=movingstd2(power,20);
