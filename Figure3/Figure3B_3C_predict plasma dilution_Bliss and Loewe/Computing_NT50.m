clc
clear all
close all

% % NT50_comp contains all the estimated NT50 values

% % Import dilution curves
Data = importdata('Fig_CI_30.mat');
Dilution(:,1) = 10.^[0:0.01:5];
Data1 = Data.epi_CI_s;

% % NT50s of all the dilution curves are compiled in 'NT50_comp' 
for i0 = 1:1:100

[~,idx]=min(abs(Data1(:,i0)-50))

if Data1(idx,i0)>=50
NT50_comp(i0,1) = interp1(Data1(idx-1:idx,i0),Dilution(idx-1:idx,1),50);
else
    Data1(idx:idx+1,i0)
    Dilution(idx:idx+1,1)
NT50_comp(i0,1) = interp1(Data1(idx:idx+1,i0),Dilution(idx:idx+1,1),50);
end
clear idx

% figure
% semilogx(Dilution(:,1),Data1(:,i0),'-',NT50_1(i0,1),[50],'*',[1 100000],[50 50],'-')
% ylim([0,100])
end

