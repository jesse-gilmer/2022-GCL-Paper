%% HEADER

% B001e_Make_Figure_1_E.m
% From: An emergent temporal basis set robustly supports cerebellar time-series learning
% Authors: Jesse I. Gilmer, Michael A. Farries, Zachary Kilpatrick, Ioannis Delis, Abigail L. Person
% 2022, BiorXiv, https://www.biorxiv.org/content/10.1101/2022.01.06.475265v1
%
% This script was written by JG.
%
% INPUTS: Make_B,C,D all need to run for this script to return correctly.
% You may comment their exection out on lines 16-18.
%
% OUTPUTS: 2 figure files to ".../Figure 1/"

%%

B001b_Make_Figure_1_B
B001c_Make_Figure_1_C
B001d_Make_Figure_1_D

clear all
close all

fdir = 'Figure 1';

a = load(fdir + "/PN_tau.mat");
a = a.tau;
b = load(fdir + "/EMG_tau.mat");
b = b.tau;
c = load(fdir + "/OU_tau.mat");
c = c.tau;

h= figure()
hold on
boxplot([a b c],[a.*0+1 b.*0+2 c.*0+3])
xticks(1:3)
xlim([.5 3.5])
xticklabels(["PN","EMG","OU"])
xlabel('Input Types')
ylabel('Decay Taus')
prettify(gcf)
saveas(h,fdir + "/Figure_1_Eun.emf")
saveas(h,fdir + "/Figure_1_Eun.png")

sc = mean(c)/100;

% a = a./sc;
% b = b./sc;
% c = c./sc;
mean(a)
std(a)
mean(b)
std(b)
mean(c)
std(c)

h= figure()
hold on
boxplot([a b c],[a.*0+1 b.*0+2 c.*0+3])
xticks(1:3)
xlim([.5 3.5])
xticklabels(["PN","EMG","OU"])
xlabel('Input Types')
ylabel('Normalied Decay Taus')
prettify(gcf)
saveas(h,fdir + "/Figure_1_E.emf")
saveas(h,fdir + "/Figure_1_E.png")