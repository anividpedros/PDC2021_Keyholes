function OE1 = opik_bplane_2_oe( theta1, phi1, zeta, xi, U, phi, Longp, ap )

ct1 = cos(theta1);
st1 = sin(theta1);

cp1 = cos(phi1);
sp1 = sin(phi1);

cp  = cos(phi);
sp  = sin(phi);


aux = 1 - U.*U - 2.*U.*ct1;

% ap = 1;
a1 = ap./ aux;
e1 = U.*sqrt( (U+2.*ct1).^2 + st1.*st1.*sp1.*sp1.*aux );
cos_in1 = (1+U.*ct1) ./sqrt(1+2.*U.*ct1 + U.*U.*(1-st1.*st1.*sp1.*sp1)) ;
% cos_in1 = (1+U.*ct1) ./(1+2.*U.*ct1 + U.*U.*(1-st1.*st1.*sp1.*sp1) );
i1 = acos( cos_in1 );

%----- Analytical Continuation ------
% Case of asteroid at r=1+d
% d = zeta.*ct1.*sp1 + xi.*cp1 ;
% a1 = a1./(1-2.*a1.*d);
% e1 = e1 + d.*U.*(2.*U - U.*st1.*st1.*sp1.*sp1 + 4.*ct1);
% tani1 = tan(i1).*(1 - d./(1+U.*ct1));
% i1 = atan(tani1);

%------------------------------------

% NEW b components
Xb = zeta.*ct1.*sp1 + xi.*cp1 ;
Yb = -zeta.*st1;
Zb = zeta.*ct1.*cp1 - xi.*sp1 ;

cos_fb1 = ( a1.*(1-e1.*e1) -ap.*(1+Xb) )./( ap.*e1.*(1+Xb) ) ;

if abs(cos_fb1)>1
    OE1 = ones(1,6)*1i ;
else


fb1 = acos( cos_fb1 ).*sign( sp );

sin_wf1 = ap.*Zb.*(1+e1.*cos(fb1)) ./ ( a1.*(1-e1.*e1).*sin(i1) );

if abs(sin_wf1) > 1
    sin_wf1 = nan;
end

wf1 = asin( sin_wf1 );
if cp < 0
    wf1 = pi - wf1;
end

w1 = wf1 - fb1;

om1_sin = sin_wf1.*cos(i1);
om1_cos = 1+cos(wf1);

Om1 = wrapTo2Pi( Longp - 2.*atan2( om1_sin , om1_cos ) + sin_wf1.*Yb.*sin(i1)./Zb );

OE1 = [a1(:) e1(:) i1(:) Om1(:) w1(:) fb1];

end