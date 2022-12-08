clc;
clear;

close all;

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
filename_dat = "sstTH.mat";
dat = sstTH;

%% parameter setting in PRA analysis
dt = dat(2,1)-dat(1,1);
orbit = [405,130,99,48.9,38,28.3,22.7,21.5,18.5];
red = 1;
pad = 5000;
srstep = 0.2;
sr1 = 0.2;
sr2 = 30;
nsim = 200;
plotn = 0;
smoothwin = 0.25;
linlog = 1;
detrended = 0;
display = 0;
method_powerspec = 1;
win = 100;
method_smooth = 'lowess';
method_cutoff = 1;
cutoff = 4;
parallelcomputing = 1;

window = 40;
step = 1;

%% detrending
if detrended ==1
    data_detrend = smooth(dat(:,1),dat(:,2),win/(dat(end,1)-dat(1,1)),method_smooth);
    dat2 = [dat(:,1),dat(:,2)-data_detrend];
    hold on;
    plot(dat2(:,1),dat2(:,2));
    detrended = 0;
else
    dat2 = dat;
end

%% ePRA 
npts = fix(window/dt);
nrow = length(dat2(:,1));
m=fix((nrow-npts)/(step/dt));
sedrate = sr1:srstep:sr2;
sedrate = sedrate';
n = length(sedrate);

corrX = zeros(m,1);
epowratioi = zeros(m,n);
epowratioi3_per = zeros(m,n);
ecyclenum = zeros(m,n);

parfor i = 1:m
    data = dat2(((i-1)*(step/dt)+1):((i-1)*(step/dt)+npts),:);
    
    [prx,cyclenum,pr,pr_sig,pr_all] = PRA(data,orbit,dt,sr1,sr2,srstep,red,nsim,pad,smoothwin,linlog,detrended,plotn,display,method_powerspec,win,method_smooth,method_cutoff,cutoff,parallelcomputing);

    corrX(i,1) = dat2(fix(npts/2+(i-1)*(step/dt)),1);
    ecyclenum(i,:) = cyclenum';
    epowratioi(i,:) = pr';
    epowratioi3_per(i,:) = pr_sig';
end


%% plot
figure
set(gcf,'unit','centimeters','position',[10,5,18,10])
set(gcf,'color','w');
ax1 = subplot(1,3,1);
[C,h]=contour(sedrate,corrX,epowratioi);
h.Fill = 'on';
u1 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u1.Label.String = 'Power ratio accumulation';
colormap(jet)
shading interp
ylabel('Depth (m)','FontSize',8,'FontName','Times New Roman')
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('powratioi','FontSize',8,'FontName','Times New Roman')
figurename ='Power ratio accumulation';
title(figurename)
set(ax1,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

ax2 = subplot(1,3,2);
[C,h]=contour(sedrate,corrX,100*(1-epowratioi3_per));
h.ShowText = 'Off';
h.Fill = 'on';
h.LevelListMode = 'manual';
colormap(jet)
shading interp
u2 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u2.Label.String = 'H_0 significance level';
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('H_0 significance level (%)','FontSize',8,'FontName','Times New Roman')
figurename ='H_0 significance level';
title(figurename)
set(ax2,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

ax3 = subplot(1,3,3);
zlevs = 0:1:length(orbit);
[C,h]=contour(sedrate,corrX,ecyclenum,zlevs);
h.Fill = 'on';
colormap(jet)
shading interp
u3 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u3.Label.String = 'Number of parameters';
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('#')
figurename = 'Number of parameters';
title(figurename)
set(ax3,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

save(filename_dat);
