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

%% data input
dat = [sstTH,sstK(:,2),sstU(:,2),sstGR(:,2),sstDT24(:,2),sstRD(:,2),sstZDEN(:,2),sstTOC(:,2)];
dt = sstTH(2,1)-sstTH(1,1);
name = ["TH","K","U","GR","DT24","RD","ZDEN","TOC"];

orbit = [405,130,99,48.9,38,28.3,22.7,21.5,18.5];
red = 1;
pad = 5000; 
srstep = 0.2;
sr1 = 0.2;
sr2 = 30;
nsim = 2000;
plotn = 1;
smoothwin = 0.25;
linlog = 1;
detrended = 0;
display = 1;
method_powerspec = 1;
win = 100;
method_smooth = 'lowess';
method_cutoff = 2;

num = length(orbit);
cutoff = rand(nsim,num)*(20-2)+2;

%% confidence levels of sensitivity ranking
[m,n] = size(dat);

data = [];

for i = 1:n-1
    if detrended == 1
        depth = dat(:,1);
        value = dat(:,i+1);
        data_detrend = smooth(depth,value,win/(depth(end)-depth(1)),method_smooth);
        dat2 = [depth,value-data_detrend];
    else
        dat2(:,1) = dat(:,1);
        dat2(:,2) = dat(:,i+1);
    end
    
    if method_powerspec == 1
        [p,f] = periodogram(dat2(:,2),[],pad,1/dt);
    else
        nw = 2;
        [p,f]=pmtm(dat2(:,2),nw,pad,1/dt);
    end
    
    P2 = p;
    
    [p,theored] = AR1noise(red,f,p,dt,smoothwin,linlog);
    data = [data,f,p];
    
end

pravalue = [];
pra_I = [];

for j = 1:nsim
    
    for i = 1:n-1
        
        [prx,nmi,pr,pr_all] = cyclepr(data(:,(i*2-1):(i*2)),orbit,sr1,sr2,srstep,method_cutoff,cutoff(j,:));
        
        pravalue(i,:) = pr;
        
    end
    
    [B,I] = sort(pravalue',2);
    pra_I(:,:,j) = I;
    
end


%% plot
sed = sr1:srstep:sr2;

RGBcolor = [255,255,255; ...
    230,245,230;...
    221,234,224; ...
    201,227,209; ...
    176,219,188;...
    126,201,146;...
    67,180,100;...
    34,139,34]/255;

figure;
set(gcf,'unit','centimeters','position',[5,2,10,10])
for z = 1:n-1
    ax = subplot(4,2,z);
    plotsum = ones(length(sed),1);
    for j = 1:n-1
        X = [sed';fliplr(sed)'];
        Y = [plotsum;zeros(length(sed),1)];
        fill(X,Y,RGBcolor(j,:),'LineStyle','none');
        ylim(ax,[0,1]);
        hold on;
        for i = 1:length(sed)
            ratio_dat(i,j) = sum(pra_I(i,j,:)==z)/nsim;
        end
        plotsum = plotsum-ratio_dat(:,j);
    end
    title(ax,name(z));
end





