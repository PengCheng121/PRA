function [prx,nmi,powratioi,powratioalli] = cyclepr(data,orbit,sr1,sr2,srstep,method_cutoff,cutoff)

sed_x = sr1:srstep:sr2;
mpts = length(sed_x);

nmi = zeros(mpts,1);
powratioi = zeros(mpts,1);
powratioalli = zeros(mpts,length(orbit));
j=1;
for i = sr1:srstep:sr2
    [nm,powratio,powratioall] = ratiovalue(data(:,1),data(:,2),orbit,i,method_cutoff,cutoff);
    nmi(j) = nm;
    powratioi(j) = powratio;
    powratioalli(j,:) = powratioall';
    j=j+1;
end

prx = linspace(sr1,sr2,j-1);
prx = prx';