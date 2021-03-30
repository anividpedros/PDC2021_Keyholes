function M = TA_2_MA( f, e)

if e < 1
    s = sqrt(1-e*e)*sin(f);
    sE = s/(1+e*cos(f));
    E = atan2( s, e+cos(f) );
    
    M = E - e*sE;
    
elseif e == 1
    error('ERR IN TA_2_MA: Parabolic case not implemented')
else
    s = sqrt(e*e-1)*sin(f);
    sH = s/(1+e*cos(f));
    H = asinh( sH );
    
    M = e*sH - H;
    
end