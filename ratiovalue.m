function [nm,powratio2,powratio] = ratiovalue(f,data_power,orbit,sr,method_cutoff,cutoff)

nm = 0;
norbits = length(orbit);

powratio = zeros(norbits,1);
powall = sum(data_power);

nyquist=f(end);
nfpts = length(f);

if method_cutoff == 1
    orbit_sed_p = 1./(orbit.*sr./100);
    orbit_max_min = zeros(2,norbits);
    
    if (orbit_sed_p(1)/cutoff) >= ((orbit_sed_p(2)-orbit_sed_p(1))/cutoff)
        orbit_max_min(1,1) = orbit_sed_p(1) - (orbit_sed_p(2)-orbit_sed_p(1))/cutoff;
        orbit_max_min(2,1) = orbit_sed_p(1) + (orbit_sed_p(2)-orbit_sed_p(1))/cutoff;
    else
        orbit_max_min(1,1) = orbit_sed_p(1) - orbit_sed_p(1)/cutoff;
        orbit_max_min(2,1) = orbit_sed_p(1) + orbit_sed_p(1)/cutoff;
    end
    
    for i = 2:norbits-1
        left_cutoff = (orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff;
        right_cutoff = (orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff;
        if(left_cutoff > right_cutoff)
            orbit_max_min(1,i) = orbit_sed_p(i)-(orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff;
            orbit_max_min(2,i) = orbit_sed_p(i)+(orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff;
        else
            orbit_max_min(1,i) = orbit_sed_p(i)-(orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff;
            orbit_max_min(2,i) = orbit_sed_p(i)+(orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff;
        end
    end
    
    orbit_max_min(1,end) =  orbit_sed_p(end)-(orbit_sed_p(end)-orbit_sed_p(end-1))/cutoff;
    orbit_max_min(2,end) =  orbit_sed_p(end)+(orbit_sed_p(end)-orbit_sed_p(end-1))/cutoff;
    
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
    
elseif method_cutoff == 2
    
    orbit_sed_p = 1./(orbit.*sr./100);
    orbit_max_min = zeros(2,norbits);
    
    if (orbit_sed_p(1)/cutoff(1)) >= ((orbit_sed_p(2)-orbit_sed_p(1))/cutoff(1))
        orbit_max_min(1,1) = orbit_sed_p(1) - (orbit_sed_p(2)-orbit_sed_p(1))/cutoff(1);
        orbit_max_min(2,1) = orbit_sed_p(1) + (orbit_sed_p(2)-orbit_sed_p(1))/cutoff(1);
    else
        orbit_max_min(1,1) = orbit_sed_p(1) - orbit_sed_p(1)/cutoff(1);
        orbit_max_min(2,1) = orbit_sed_p(1) + orbit_sed_p(1)/cutoff(1);
    end
    
    for i = 2:norbits-1
        left_cutoff = (orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff(i);
        right_cutoff = (orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff(i);
        if(left_cutoff > right_cutoff)
            orbit_max_min(1,i) = orbit_sed_p(i)-(orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff(i);
            orbit_max_min(2,i) = orbit_sed_p(i)+(orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff(i);
        else
            orbit_max_min(1,i) = orbit_sed_p(i)-(orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff(i);
            orbit_max_min(2,i) = orbit_sed_p(i)+(orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff(i);
        end
    end
    
    orbit_max_min(1,end) =  orbit_sed_p(end)-(orbit_sed_p(end)-orbit_sed_p(end-1))/cutoff(end);
    orbit_max_min(2,end) =  orbit_sed_p(end)+(orbit_sed_p(end)-orbit_sed_p(end-1))/cutoff(end);
    
else

    orbit_max_min = 1./(cutoff.*sr./100);

end

for i = 1 : norbits
    
    if nyquist <= orbit_sed_p(i) || f(2) >= orbit_sed_p(i)
        nm = nm + 1;
        continue;
    end
    
    nfmin=ceil(nfpts*orbit_max_min(1,i)/nyquist);
    
    nfmax=fix(nfpts*orbit_max_min(2,i)/nyquist);
    
    if orbit_max_min(2,i)/nyquist>1
        nfmax = nfpts;
    end
    
    nfn=nfmax-nfmin+1;
    spq=zeros(1,nfn);
    ij=1;
    for q=nfmin:nfmax
        spq(ij)=data_power(q);
        ij=ij+1;
    end
    
    powratio(i,1) = sum(spq)/powall;
    
    if sum(spq) <= 0
        nm = nm + 1;
    end
    
end

powratio2 = sum(powratio);










