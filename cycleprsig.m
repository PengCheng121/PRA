function [powratioi,powratioalli] = cycleprsig(data,orbit,sr1,sr2,srstep,method_cutoff,cutoff)

sed_x = sr1:srstep:sr2;  % tested sed. rate series
mpts = length(sed_x);  % tested sed. rates number

powratioi = zeros(mpts,1);
powratioalli = zeros(mpts,length(orbit));
j=1;
for i = sr1:srstep:sr2
    [~,powratio,powratioall] = ratiovalue(data(:,1),data(:,2),orbit,i,method_cutoff,cutoff);
    powratioi(j) = powratio;
    powratioalli(j,:) = powratioall';
    j=j+1;
end
