% Keplerian elements into structure for MOID fxn
function A = K2S(OE,AU) 
A.sma   = OE(1)/AU;
A.e     = OE(2);
A.i     = OE(3);
A.Omega = OE(4);
A.argp  = OE(5);
end