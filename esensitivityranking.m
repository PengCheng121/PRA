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
nsim = 100;
plotn = 1;
smoothwin = 0.25;
linlog = 1;
detrended = 0;
display = 1;
method_powerspec = 1;
win = 100;
method_smooth = 'lowess';
method_cutoff = 2;

window = 40;
step = 1;

%% confidence levels of sensitivity ranking
npts = fix(window/dt);
nrow = length(dat(:,1));
mm=fix((nrow-npts)/(step/dt));
sedrate = sr1:srstep:sr2;
sedrate = sedrate';
n = length(sedrate);

corrX = zeros(mm,1);
epowratioi = zeros(mm,n);

hwaitbar = waitbar(0,'Processing ... [CTRL + C to quit]','WindowStyle','normal','Name','Wait Bar');

for z = 1:mm
    data = dat(((z-1)*(step/dt)+1):((z-1)*(step/dt)+npts),:);
    
    num = length(orbit);
    
    cutoff = rand(nsim,num)*(20-2)+2;
    
    
    
    [m,n] = size(data);
    
    data2 = [];
    
    for i = 1:n-1
        if detrended == 1
            depth = data(:,1);
            value = data(:,i+1);
            data_detrend = smooth(depth,value,win/(depth(end)-depth(1)),method_smooth);
            dat2 = [depth,value-data_detrend];
        else
            dat2(:,1) = data(:,1);
            dat2(:,2) = data(:,i+1);
        end
        
        if method_powerspec == 1
            [p,f] = periodogram(dat2(:,2),[],pad,1/dt);
        else
            nw = 2;
            [p,f]=pmtm(dat2(:,2),nw,pad,1/dt);
        end
        
        P2 = p;
        
        %% remove AR1 noise
        [p,theored] = AR1noise(red,f,p,dt,smoothwin,linlog);
        data2 = [data2,f,p];
        
    end
    
    pravalue = [];
    pra_I = [];
    for j = 1:nsim
        
        for i = 1:n-1
            
            %% MPDA
            
            [prx,nmi,pr,pr_all] = cyclepr(data2(:,(i*2-1):(i*2)),orbit,sr1,sr2,srstep,method_cutoff,cutoff(j,:));
            
            pravalue(i,:) = pr;
            
        end
        
        [B,I] = sort(pravalue',2);
        pra_I(:,:,j) = I;
        
    end
        
    
    
    
    corrX(z,1) = dat(fix(npts/2+(z-1)*(step/dt)),1);
    
    for j = 1:length(sedrate)
        ratio(1,j) = sum(pra_I(j,end,:)==1)/nsim;
    end
    epowratioi(z,:) = ratio;
    
    disp(['>> Simulation ',num2str(z),' of ',num2str(mm)]);
    waitbar(z/mm);
    
end

if ishandle(hwaitbar)
    close(hwaitbar);
end

%% plot
figure
set(gcf,'unit','centimeters','position',[10,5,18,10])
set(gcf,'color','w');
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
set(gca,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');















