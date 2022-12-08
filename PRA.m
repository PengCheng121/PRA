function [prx,cyclenum,pr,pr_sig,pr_all] = PRA(dat,orbit,dt,sr1,sr2,srstep,red,nsim,pad,smoothwin,linlog,detrended,plotn,display,method_powerspec,win,method_smooth,method_cutoff,cutoff,parallelcomputing) 
% power ratio accumulation, or PRA
%
% input parameters
%     dat:  Two columns of data. The first colume is depth. The second colume is depth (Unit is m). The second column is value.
%     orbit: target orbital parameters
%     dt: sampling rate of data (dat)
%     sr1: begining sedimentation rates to be estimated (Unit is cm/kyr).
%     sr2: end sedimentation rates to be estimated (Unit is cm/kyr).
%     srstep: step of sedimentation rates (Unit is cm/kyr).
%     red: 0 = no remove red noise. 1 = remove red noise.
%     nsim: number of Monte Carlo simulation
%     pad: zero-padding
%     smoothwin: smoothing window
%     linlog: fit to S(f) or logS(f). 1 = linear; 2 = log
%     detrended: 1 = detrending. else = no detrending
%     plotn: 1 = plot power spectra. else = no plot
%     display: 1 = show wait bar. else = no show.
%     method_powerspec: 1 = Periodogram power spectral density estimate. 2 = Multitaper power spectral density estimate.
%     win: Window length in detrending
%     method_smooth: method of smoothing
%     method_cutoff: 1 = "cutoff" is number of intervals. else(except 2) = "cutoff" is cutoff frequencies.
%     cutoff: cutoff depends on method_cutoff.
%     parallelcomputing: 0 = no parallel for loop. 1 = parallel for loop.
%
% output parameters
%     prx: sedimentation rate
%     cyclenum: Number of parameters
%     pr: Cumulative power ratio
%     pr_sig :H0 significant level
%     pr_all: power ratio

if nargin < 20
    parallelcomputing = 0;
    if nargin < 19
        cutoff = 4;
        if nargin < 18
            method_cutoff = 1;
            if nargin < 17
                method_smooth = 'lowess';
                if nargin < 16
                    win = (dat(end,1)-dat(1,1))*0.35;
                    if nargin < 15
                        method_powerspec = 1;
                        if nargin < 14
                            display = 1;
                            if nargin < 13
                                plotn = 1;
                                if nargin < 12
                                    detrended = 1;
                                    if nargin < 11
                                        linlog = 1;
                                        if nargin < 10
                                            smoothwin = 0.25;
                                            if nargin < 9
                                                pad = 5000;
                                                if nargin < 8
                                                    nsim = 1000;
                                                    if nargin < 7
                                                        red = 3;
                                                        if nargin < 6
                                                            error('Too few input arguments');
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% detrending
if detrended == 1
    depth = dat(:,1);
    value = dat(:,2);
    data_detrend = smooth(depth,value,win/(depth(end)-depth(1)),method_smooth);
    dat2 = [depth,value-data_detrend];
else
    dat2 = dat;
end

%% power spectral analysis
if method_powerspec == 1
    [p,f] = periodogram(dat2(:,2),[],pad,1/dt);
else
    nw = 2;
    [p,f]=pmtm(dat2(:,2),nw,pad,1/dt);
end

P2 = p;

%% remove red noise
[p,theored] = AR1noise(red,f,p,dt,smoothwin,linlog);
data = [f,p];

%% plot power spectra
if plotn == 1
    
    figure;
    AX1 = subplot(3,1,1);
    plot(AX1,dat(:,1),dat(:,2),'b','LineWidth',1);
    hold on;
    if detrended == 1
        plot(AX1,dat2(:,1),dat2(:,2),'r','LineWidth',1);
        legend(AX1,'Data series','Detrended data series');
    else
        legend(AX1,'Data series');
    end
    title(AX1,'Data series');
    ylabel(AX1,'Value');
    xlabel(AX1,'Depth (m)');
    set(AX1,'XMinorTick','on','YMinorTick','on');
    
    AX2 = subplot(3,1,2);
    plot(AX2,f,P2,'b','LineWidth',1);
    hold on;
    plot(AX2,f,theored,'r','LineWidth',1);
    title(AX2,'Raw periodogram and red noise');
    xlabel(AX2,'Frequency (cycle/m)');
    ylabel(AX2,'Power');
    legend(AX2,'Power spectrum of data series','AR(1) noise');
    set(AX2,'XMinorTick','on','YMinorTick','on');
    
    AX3 = subplot(3,1,3);
    plot(AX3,f,p,'b','LineWidth',1);
    xlabel(AX3,'Frequency (cycle/m)');
    ylabel(AX3,'Power');
    title(AX3,'Red noise removed');
    legend(AX3,'Power spectrum of data series');
    set(AX3,'XMinorTick','on','YMinorTick','on');

end

%% PRA
[prx,nmi,pr,pr_all] = cyclepr(data,orbit,sr1,sr2,srstep,method_cutoff,cutoff);
mpts = length(prx);
orbitn = length(orbit);
cyclenum = (orbitn-nmi);

%% Monte Carlo simulation
n = length(dat2(:,1));
rho = rhoAR1(dat2(:,2));
yred=filter(1,[1;-rho],randn(n,nsim));
pr_mc = zeros(mpts,nsim);

if parallelcomputing == 0
    
    if display == 1
        hwaitbar = waitbar(0,'Monte Carlo processing ... [CTRL + C to quit]','WindowStyle','normal','Name','Wait Bar');
    end
    
    for i = 1: nsim
        if method_powerspec == 1
            [randspectrum,f] = periodogram(yred(:,i),[],pad,1/dt);
        else
            nw = 2;
            [randspectrum,f]=pmtm(yred(:,i),nw,pad,1/dt);
        end
        
        [randspectrum_p,~] = AR1noise(red,f,randspectrum,dt,smoothwin,linlog);
        sim_spectum = [f,randspectrum_p];
        
        [pr_rand,~] = cycleprsig(sim_spectum,orbit,sr1,sr2,srstep,method_cutoff,cutoff);
        if display == 1
            disp(['>> Simulation ',num2str(i),' of ',num2str(nsim)]);
            waitbar(i/nsim);
        end
        pr_mc(:,i) = pr_rand;
    end
    
    if display == 1
        if ishandle(hwaitbar)
            close(hwaitbar);
        end
    end
    
else
    
    parfor i = 1: nsim
        
        if method_powerspec == 1
            [randspectrum,f] = periodogram(yred(:,i),[],pad,1/dt);
        else
            nw = 2;
            [randspectrum,f]=pmtm(yred(:,i),nw,pad,1/dt);
        end
        
        [randspectrum_p,~] = AR1noise(red,f,randspectrum,dt,smoothwin,linlog);
        sim_spectum = [f,randspectrum_p];
        
        [pr_rand,~] = cycleprsig(sim_spectum,orbit,sr1,sr2,srstep,method_cutoff,cutoff);

        pr_mc(:,i) = pr_rand;
    end
    
end

%% MC results
pr_sig = zeros(mpts,1);

for i = 1: mpts
    pr_sim = pr_mc(i,:);
    pr2 = pr(i);
    
    pr_sim2 = pr_sim(pr_sim<pr2);
    pr_totallength = length(pr_sim(~isnan(pr_sim)));
    
    pr_sig(i) = (pr_totallength-length(pr_sim2)+1)/(pr_totallength+1);
    if pr_sig(i) == 0
        pr_sig(i) = 1/(pr_totallength+1);
    end
end



