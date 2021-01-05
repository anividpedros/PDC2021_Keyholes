function Eq = Kep_2_Equinoctial( Kep )

e  = Kep(2) ;
ti = Kep(3) ;
w  = Kep(5)+Kep(4) ;
th = Kep(4) ;

h = e*sin( w );
k = e*cos( w );
p = ti*sin( th );
q = ti*cos( th );

Eq = [ Kep(1); h; k; p; q; Kep(6) ];