function [p,theored] = AR1noise(red,f,p,dt,smoothwin,linlog)

if red == 0
    theored = zeros(length(p),1);

elseif red == 1
    % robust
    theored = redconf_any(f,p,dt,smoothwin,linlog);
    p = p - theored;
    p(p<0) = 0;   

end