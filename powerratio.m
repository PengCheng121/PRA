clc;
clear;

close all;

tic;

filepath = "D:\matlab2019a\bin\work\PRAmethod\Data\";

%% load data
filename = filepath + "sk1sK-100-LOWESS.txt";
load(filename)
sstK = sk1sK_100_LOWESS;

filename = filepath + "sk1sGR-100-LOWESS.txt";
load(filename)
sstGR = sk1sGR_100_LOWESS;

filename = filepath + "sk1sTH-100-LOWESS.txt";
load(filename)
sstTH = sk1sTH_100_LOWESS;

filename = filepath + "sk1sU-100-LOWESS.txt";
load(filename)
sstU = sk1sU_100_LOWESS;

filename = filepath + "sk1sDT24-100-LOWESS.txt";
load(filename)
sstDT24 = sk1sDT24_100_LOWESS;

filename = filepath + "sk1slog10RD-100-LOWESS.txt";
load(filename)
sstRD = sk1slog10RD_100_LOWESS;

filename = filepath + "sk1sZDEN_1025.13_1128.17-100-LOWESS.csv";
load(filename)
sstZDEN = sk1sZDEN_1025_13_1128_17_100_LOWESS;

filename = filepath + "sk1sknlogRtoc_1025.13_1128.17-100-LOWESS.csv";
load(filename)
sstTOC = sk1sknlogRtoc_1025_13_1128_17_100_LOWESS;


%% data input in PRA analysis
filename_dat = "sstTH.csv";
dat = sstTH;

%% parameter setting in PRA analysis
dt = dat(2,1)-dat(1,1);
orbit = [405,130,99,48.9,38,28.3,22.7,21.5,18.5];
red = 1;
pad = 5000;
srstep = 0.2;
sr1 = 0.2;
sr2 = 30;
nsim = 500;
plotn = 1;
smoothwin = 0.25;
linlog = 1;
detrended = 0;
display = 1;
method_powerspec = 1;
win = 100;
method_smooth = 'lowess';
method_cutoff = 1;
cutoff = 4;
parallelcomputing = 1;

%% PRA analysis
[prx,cyclenum,pr,pr_sig,pr_all] = PRA(dat,orbit,dt,sr1,sr2,srstep,red,nsim,pad,smoothwin,linlog,detrended,plotn,display,method_powerspec,win,method_smooth,method_cutoff,cutoff,parallelcomputing);

%% save results
Header = {'Sedimentation rate','Number of parameters','Cumulative power ratio','H0 significant level'};
n = length(orbit);
for i = 1:n
    Header{end+1} = num2str(orbit(i));
end
result1 = num2cell([prx,cyclenum,pr,pr_sig,pr_all]);
outputdata = [Header;result1]; 
xlswrite(filename_dat,outputdata);

%% plot
RGBcolor = [254,67,101;...
    252,157,154;...
    249,205,173;...
    217,116,43;...
    210,180,140;...
    230,180,80;...
    61,145,64;...
    137,157,192;...
    174,221,129;...
    131,175,155;...
    127,255,212;...
    153,77,82;...
    200,200,169]/255;

mpts = length(prx);
figure;
set(gcf,'unit','centimeters','position',[10,5,7.5,7.5])
set(gcf,'color','w');
ax1 = subplot('Position',[0.15 0.7 0.55 0.25]);
plotsum = pr;
for i = length(orbit):-1:1
    X = [prx;fliplr(prx')'];
    Y = [plotsum;zeros(mpts,1)];
    fill(ax1,X,Y,RGBcolor(i,:),'LineStyle','none');
    hold on;
    plotsum = plotsum-pr_all(:,length(orbit)-i+1);
    leg{i} = num2str(orbit(i));
end
lgd = legend(leg,'NumColumns',1,'FontSize',8,'Location','NorthEastOutside','FontName','Times New Roman','Position',[0.85 0.7 0 0]);
ylabel(ax1,'Power ratio','FontSize',8,'FontName','Times New Roman');
set(ax1,'XMinorTick','on','FontSize',8,'FontName','Times New Roman','xticklabel',[]);

ax2 = subplot('Position',[0.15 0.4 0.55 0.25]);
semilogy(ax2,prx,pr_sig,'r','LineWidth',1);
ylabel(ax2,'H_0 significance level','FontSize',8,'FontName','Times New Roman');
ylim(ax2,[0.5*min(pr_sig) 1]);
line([sr1, sr2],[.10, .10],'LineStyle',':','Color','k');
line([sr1, sr2],[.05, .05],'LineStyle',':','Color','k');
line([sr1, sr2],[.01, .01],'LineStyle','--','Color','k');
line([sr1, sr2],[.001, .001],'LineStyle',':','Color','k');
set(ax2,'XMinorTick','on','FontSize',8,'FontName','Times New Roman','xticklabel',[],'Ydir','reverse');

ax3 = subplot('Position',[0.15 0.1 0.55 0.25]);
plot(ax3,prx,cyclenum,'b','LineWidth',1);
xlabel(ax3,'Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman');
ylabel(ax3,{'Number of', 'parameters'},'FontSize',8,'FontName','Times New Roman');
ylim(ax3,[0 length(orbit)+0.5]);
set(ax3,'XMinorTick','on','FontSize',8,'FontName','Times New Roman');

a = sgtitle(filename_dat,'FontSize',8,'FontName','Times New Roman');


toc;





