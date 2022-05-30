function [co2 Hnew pH k1 k2 ff]= co2_follows(T,S,dic,pt,sit,ta,pHlocal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE co2 and new pH for given DIC.
% Efficient solver following Follows et al (2006), 
% On the solution of the Carbonate chemistry system
% in ocean biogeochemistry models, Ocean modelling,12,290-301
%       co2 = dissolved inorganic carbon (mol./kg)
%       dic = total inorganic carbon (mol./kg) 
%       ta  = total alkalinity (mol./kg) 
%       pt  = inorganic phosphate (mol./kg) 
%       sit = inorganic silicate (mol./kg) 
%       T   = temperature (degrees K)
%       S   = salinity (PSU)
%		hg  = first guess of [H+] 

% set first guess for [H]
    hg = 10.^(-pHlocal);
% estimated concentration of borate based on salinity
    scl=S./1.80655;
    bt=0.000232 .* scl./10.811;
% some definitions ...
    S2 = S.*S ;
    sqrtS = S.^0.5 ;
    invT = 1.0./T ;
    T1 = T./100.0 ;

% Coefficient algorithms as used in OCMIP2 protocols 
% K1, K2 Millero (1995) using Mehrbach data
k1 = 10.^(-1.*( 3670.7.*invT - ...
	62.008 + 9.7944.*log(T) - ...
	0.0118.*S + 0.000116.*S2)) ;

k2 = 10.^(-1.*(1394.7.*invT + 4.777 - ...
	0.0184.*S + 0.000118.*S2)) ;

% K1p, K2p, K3p  DOE (1994) 
k1p = exp(-4576.752.*invT + 115.525 - ...
	18.453.*log(T) + ...
	(-106.736.*invT + 0.69171).*sqrtS + ...
	(-0.65643.*invT - 0.01844).*S) ;

k2p = exp(-8814.715.*invT + 172.0883 - ...
	27.927.*log(T) + ...
	(-160.34.*invT + 1.3566).*sqrtS + ...
	(0.37335.*invT - 0.05778).*S) ;

k3p = exp(-3070.75.*invT - 18.141 + ...
	(17.27039.*invT + 2.81197) .* ...
	sqrtS + (-44.99486.*invT - 0.09984).*S) ;

% Kb, Millero (1995) using data from Dickson
kb = exp((-8966.90 - 2890.53.*sqrtS - 77.942.*S + ...
	1.728.*S.^1.5 - 0.0996.*S2).*invT + ...
	(148.0248 + 137.1942.*sqrtS + 1.62142.*S) + ...
	(-24.4344 - 25.085.*sqrtS - 0.2474.*S) .* ...
	log(T) + 0.053105.*T.*sqrtS) ;

% Kw, Millero (1995)
kw = exp(-13847.26.*invT + 148.9652 - ...
	23.6521.*log(T) + ...
	(118.67.*invT - 5.977 + 1.0495.*log(T)) .* ...
	sqrtS - 0.01615.*S) ;

% Ksi, Millero (1995)
I = (19.924.*S)./(1000 - 1.005.*S) ;
ksi = exp(-8904.2.*invT + 117.385 - ...
	19.334.*log(T) + ...
       	(-458.79.*invT + 3.5913).*(I.^.5) + ...
       	(188.74.*invT - 1.5998).*I + ...
	(-12.1652.*invT + 0.07871).*(I.*I) + ...
	log(1.0 - 0.001005.*S)) ;

% fugacity, Weiss and Price, Marine Chem, 8, 347 (1990)
 ff = exp(-162.8301 + 218.2968./(T1) + ...
	90.9241.*log(T1) - 1.47696.*(T1.*T1) + ...
	S.*(.025695 - .025225.*(T1) + ...
	0.0049867.*(T1.*T1)));

% First guess of [H+]: from last timestep .*OR.* fixed for cold start
% --- here iterate for accurate solution
for i = 1:20
	% estimate contributions to total alk from borate, silicate, phosphate
	bohg = (bt.*kb)./(hg+kb);
	siooh3g = (sit.*ksi)./(ksi + hg);
	denom = (hg.*hg.*hg) + (k1p.*hg.*hg) + (k1p.*k2p.*hg) + (k1p.*k2p.*k3p);
	h3po4g = (pt.*hg.*hg.*hg)./denom;
	h2po4g = (pt.*k1p.*hg.*hg)./denom;
	hpo4g = (pt.*k1p.*k2p.*hg)./denom;
	po4g = (pt.*k1p.*k2p.*k3p)./denom;
	% estimate carbonate alkalinity
	cag = ta - bohg - (kw./hg) + hg - hpo4g - 2.*po4g + h3po4g - siooh3g ;
	% estimate hydrogen ion conc
        gamm=dic./cag;
        dummy=(1-gamm).^2.*k1.^2-4.*k1.*k2.*(1-2.*gamm);
        Hnew = 0.5.*((gamm-1).*k1+sqrt(dummy)); 
        hg = Hnew ;
end
% evaluate co2
 co2 = dic./(1.0 + (k1./Hnew) + (k1.*k2./(Hnew.*Hnew))) ;

% calc final pH
pH = -log10(Hnew) ;
end

