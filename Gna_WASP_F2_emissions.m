%Katherine Turner 2022 Earth's Future
%In this version the atmosphere-ocean model contains an additional
%terrestrial model in the style of WASP (Goodwin et al., 2016) that
%contains two boxes for soil and vegetation carbon and linearised
%sensitivities to temperature and atmospheric CO2 changes

%Anna Katavouta 2019 Journal of Climate (based on Gnanadesikan 1999)
%!!!!!IN THIS VERSION F2: For carbon flux between atmosphere-ocean the
%atmosphere sees the ocean surface as one region and in this case the DIC
%and temperature of the ocean that the atmosphere sees is a WEIGHTED AREA AVERAGE of 
%the boxes in contact with the atmosphere (i.i, mixed layer in
%tropics-subtropics, upper layer in northern high latitudes and southern
%ocean).
% THE ANTHROPOGENIC HEAT AND CARBON flux from the atmosphere into the ocean is
% distributed UNEVENLY into the ocean (for example higher into the Southern Ocean, choose as you want:
% Fa_s Fa_n Fa_m Hm Hs Hn in lines 251-254 and 283)

tau=0.08;%wind stress in southern ocean (N/m2)
FRAC=0.5;% delta fraction of isolation: roughly capture the percentage of subduction 
         % occuring in the Southern Ocean (values between 0 and 1)
lambda=0.6;%climate feedback
dNPP_dT=-2;% temperature dependence of NPP, range -5 to 1 PgC yr-1 K-1
gamma_co2=0.5;% co2 fertilisation rate, range 0-1
dtau_dT=-0.5;% dependency of soil residency time on temperature, range -2 to 1 yr K-1

RCP = load('/Users/keturner/Gnanadesikan_model_comments/global_alkalinity/G2_fcts/forcings.mat', 'Iem');
% choose RCP here- all run from 1850 to 2300 (450 years)
% 1 = RCP 4.5
% 2 = RCP 6.0
% 3 = RCP 8.5
% 4 = RCP 2.6
Iem=RCP.Iem(4).data(1:end-1);

Rate_weath=2;%weathering rate in Pg C/year
halftime = 1; %set to 1 if you want weathering to occur at half the rate but for twice the time

time = 1850 + (0:164249)/365;

Rate_weath_ext = zeros(1, 450*365);
Iweath = zeros(1, 450*365);
Rate_weath_ext(1:150*365) = 0;
Iweath(1:150*365) = 0;
if halftime==0
    Rate_weath_ext(150*365+1:250*365)=Rate_weath./365;
    Rate_weath_ext(250*365+1:450*365)=0;
    Iweath(150*365+1:250*365)=Rate_weath.*((1:100*365)./365);
    Iweath(250*365+1:450*365)=Iweath(250*365);
end
if halftime==1
    Rate_weath_ext(150*365+1:350*365)=Rate_weath./365.*0.5;
    Rate_weath_ext(350*365+1:450*365)=0;
    Iweath(150*365+1:350*365)=Rate_weath.*0.5.*((1:200*365)./365);
    Iweath(350*365+1:450*365)=Iweath(350*365);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% general parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=5.1352*10^18;%mass of dry air
Kg=5*10^(-5); %air-sea gas transfer coefficient
P=1; %atmospheric pressure in atm
Ma=28.97*10^(-3);% mean molecular mass
Na=m./Ma; %number of moles in the atmosphere
alpha=5.35;%radiative forcing coefficient
CC=50;%air-sea heat transfer parameter 

rho_o=1025;%referenced ocean density
rho_a=1;% referenced air density
cp_o=4*10^3;%ocean specific heat capacity
cp_a=1*10^3;% air specific heat capacity
Sref=34.5;%referenced ocean salinity

h_a=10*10^3;%thickness of atmospheric layer
Ar_a=3.4*10^14;%area of atmospheric layer
Ar_o=Ar_a;% area of ocean
Ar_gna=Ar_a*0.6;%area of tropics+subtropic
Ar_s=(Ar_o-Ar_gna)./2;% Area of Southern ocean
Ar_n=Ar_s;% Area of Northern high latitudes

D=4000;%depth of the whole ocean
h_m=100;%depth of mixed layer 
D_n=1000;% depth of the Northern high latitudes upper box (convection depth-like)
D_s=1000;% depth of the Southern Ocean upper box

thetam=25+273.15;%reference temperature of mixed layer
thetad=5+273.15;%reference temperature of deep box
thetas=5+273.15;%reference temperature of southern ocean upper box
thetan=5+273.15;%reference temperature of northern high latitudes upper box

Pho=2.157*10^(-6);%mol/kg Inorganic phosphate for the ocean; choose value; will remain unchanged 
Sil=2.157*10^(-6)*32;%mol/kg Silicate value for the ocean; choose value; will remain unchanged

Deltat=24*3600; %time step: a day in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters for the terrestrial boxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NPP0 = 60; % net primary production, PgC yr−1
LL0 = 60; % leaf litter, PgC yr−1
SR0 = 60; % soil respiration, PgC yr−1
Iveg0 = 450; % PgC
Isoil0 = 1500; % PgC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters for the dynamic overturning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimate wind driven Ekman transport in the southern ocean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx=3*10^7; %(m) zonal extend of Southern ocean
fc=10^(-4);%(1/s)coriolis parameter
VTek=(tau*Lx)/(rho_o*fc);% transport (m^3/s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimate eddy driven component in the southern ocean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ked=1000;%(m^2/s) horizontal eddy diffusivity
Ly=1.5*10^6; % (m) width of ACC current
VTedh=Ked*Lx/Ly;% (m^2/s) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%North Atlantic deep water formation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=0.02;%m/s2 reduced gravity
VTnah=g/(2*fc);%(m/s) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diapycnal transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kd=1*10^(-5);%(m^2/s) diapycnal mixing coefficient
VTuh=Ar_gna*Kd;  %(m^4/s) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate upper ocean depth (thermocline+mixed layer) for a steady sate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%no transport changes in time so all the transport in and out of the box
%should be zero %Tek-Ted-Tna+Tu=0.
p=[-VTnah -VTedh VTek VTuh];
DD=roots(p);
for ii=1:length(DD)
    if DD(ii)<0
        DD(ii)=nan;
    end
    if imag(DD(ii))>0
        DD(ii)=nan;
    end
end
DD(isnan(DD))=[];
DD(DD==0)=[];
% where DD is the depth of the mixed laye+thermocline (upper ocean in
% subtropics+tropics)
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate the transport that corresponds to the upper ocean depth in steady state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%assume positive into the upper ocean and negatice into the deep ocean
VTek=VTek;% Ekman
VTu=VTuh/DD;% diapycnal
VTna=-VTnah*DD^2;% sinking northern high latitudes
VTed=-VTedh*DD;% eddy return flow in Southern Ocean
VTRO=VTek+VTed;% this is the combined southern ocean transport (Ekman+eddy return flow)

VTu+VTek+VTna+VTed;%should give zero if steady state else the changes in volume
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate heat fluxes at the pre-industrial in the steady state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_s=thetad;T_n=thetad;T_d=thetad;T_m=thetam;
T_u=-(FRAC*VTRO*T_s+(1-FRAC)*VTRO*T_m+VTu*T_d)./VTna;% estimate thermocline temperature at steady state
Hm_pre=-(rho_o.*cp_o).*((1-FRAC)*VTRO*T_s-(1-FRAC)*VTRO*T_m);% (J/s) heat into mixed layer at tropics subtropics 
Hn_pre=-(rho_o.*cp_o)*(-VTna*T_u+VTna*T_n);%(J/s) heat into the upper ocean at the northern high latitudes
%pre-industrial flux at southern ocean is assumed zero

% T_cont is the surface temperature that the atmosphere sees SST kind 
T_cont=(T_s.*Ar_s+T_n.*Ar_n+T_m.*Ar_gna)./Ar_o;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate carbon fluxes at the pre-industrial in the steady state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chi=280*10^(-6);%preindustrial mixing ratio
At=2350*10^(-6);%pre-indutrial alkalinity

% 16 July 2019 -- set up all regional alkalinities to be the same
alk_s(1)=At; alk_n(1)=At; alk_m(1)=At; alk_u(1)=At; alk_d(1)=At;

%To have steady state estimate what will be the surface DIC in equilibrium with the pco2
[dic_eq Heq pHeq K1m K2m Kom]= Equilibr(T_cont,Sref,chi,Pho,Sil,At,8);
%set the value of the high latitudes and then find the corresponding value
%at the tropics-subtropics for equilibrium: essentially the area weighted
%average should be equal to dic_eq.
dic_n_pre=2140*10^(-6);
dic_s_pre=dic_n_pre;
dic_m_pre=(dic_eq*(Ar_s+Ar_n+Ar_gna)-dic_n_pre*Ar_n-dic_s_pre*Ar_s)./(Ar_gna);

dic_d_pre=dic_n_pre;
dic_u_pre=-(FRAC*VTRO*dic_s_pre+(1-FRAC)*VTRO*dic_m_pre+VTu*dic_d_pre)./VTna;%estimate thermocline DIC at steady state
Fa_m_pre=-(rho_o./Ar_gna).*((1-FRAC)*VTRO*dic_s_pre-(1-FRAC)*VTRO*dic_m_pre);% mol/(m^2*s) carbon into mixed layer at tropics subtropics 
Fa_n_pre=-(rho_o./Ar_n)*(-VTna*dic_u_pre+VTna*dic_n_pre);% mol/(m^2*s) carbon into mixed layer at north high latitudes
%pre-industrial flux at southern ocean is assumed zero

% dic_carb here is the surface DIC that the atmosphere sees: you can choose different
% definitions but here I define it as as an area weighting of the DIC
% in the boxes in contact with the atmosphere
dic_carb_pre=(dic_s_pre.*Ar_s+dic_n_pre.*Ar_n+dic_m_pre.*Ar_gna)./Ar_o;

% define pH etc. at the pre-industrial
[co2new H1 pH1,K1,K2,Ko]= co2_follows(T_cont,Sref,dic_carb_pre,Pho,Sil,At,8);
H_carb_pre=H1;
pH_carb_pre=pH1;
CO2_carb_pre=dic_carb_pre./(1+(K1/H_carb_pre)+(K1.*K2./(H_carb_pre^2)));
Fa_test=-rho_o.*Kg.*(CO2_carb_pre-Ko.*P.*chi);% this should be close to zero 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXTRA STEP TO CHECK YOU START FROM A A STEADY STATE AT PRE-INDUSTRIAL
% Use it only if you want but not needed for the actual run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chi_pre=chi;
% %check if start from equilibrium
% for tt=2:5000*365
%     chi_pre(tt)=chi_pre(tt-1);
%     dic_s_pre(tt)=dic_s_pre(tt-1)+(Deltat./(Ar_s*D_s))*(VTRO*dic_d_pre(tt-1)-VTRO*dic_s_pre(tt-1));
%     dic_u_pre(tt)=dic_u_pre(tt-1)*((DD-h_m)./(DD-h_m))...
%         +(Deltat./(Ar_gna*(DD-h_m)))*(FRAC*VTRO*dic_s_pre(tt-1)+(1-FRAC)*VTRO*dic_m_pre(tt-1)+VTna*dic_u_pre(tt-1)+VTu*dic_d_pre(tt-1));
%     dic_n_pre(tt)=dic_n_pre(tt-1)+(Deltat./(Ar_n*D_n.*rho_o))*(Fa_n_pre(tt-1))*Ar_n+(Deltat./(Ar_n*D_n))*(-VTna*dic_u_pre(tt-1)+VTna*dic_n_pre(tt-1));
%     dic_m_pre(tt)=dic_m_pre(tt-1)+Deltat./(Ar_gna*h_m*rho_o)*(Fa_m_pre(tt-1))*Ar_gna+Deltat./(Ar_gna*h_m)*((1-FRAC)*VTRO*dic_s_pre(tt-1)-(1-FRAC)*VTRO*dic_m_pre(tt-1));
%     dic_d_pre(tt)=dic_d_pre(tt-1)*((Ar_gna*(D-DD)+Ar_s*(D-D_s)+Ar_n*(D-D_s))./(Ar_gna*(D-DD)+Ar_s*(D-D_s)+Ar_n*(D-D_s)))...
%               +(Deltat./(Ar_gna*(D-DD)+Ar_s*(D-D_s)+Ar_n*(D-D_s)))*(-VTna*dic_n_pre(tt-1)-VTRO*dic_d_pre(tt-1)-VTu*dic_d_pre(tt-1));
% 
% Fa_m_pre(tt)=-(rho_o./Ar_gna).*((1-FRAC)*VTRO*dic_s_pre(tt)-(1-FRAC)*VTRO*dic_m_pre(tt));
% Fa_n_pre(tt)=-(rho_o./Ar_n)*(-VTna*dic_u_pre(tt)+VTna*dic_n_pre(tt));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set your initial values for the simulation with emissions
%  (some of the initial values at zero at pre-industrial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dic_m=dic_m_pre(end);dic_d=dic_d_pre(end);
dic_u=dic_u_pre(end);dic_s=dic_s_pre(end);dic_n=dic_n_pre(end);
Fa_all=Fa_test;Fa_all_pre=Fa_test;pH_carb=pH_carb_pre;
Iveg=Iveg0;Isoil=Isoil0;
%Initial values for the carbon and heat fluxes associated with emissions (should be
%zero in the pre-industrial as there are no emissions)
Fa_s=0;Fa_n=0;Fa_m=0;
dI_v=0;dI_s=0;
Hn=0;Hs=0;Hm=0;
NTOA=0;N=0;DTa=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Experiment with carbon emissions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% new scenario with weathering

%16 July 2019- included weathering emissions. We go ahead and convert the
%yearly fluxes into surface fluxes that are evenly partitioned between all
%surface boxes
%THIS IS SAME AS IN GNA_F_1
Ga = Rate_weath_ext*(1/(24*60^2))*1e15/12.011/Ar_o; % weathering rate in mol C (moleq) / m^2 /s
Ga_s=Ga; Ga_n=Ga; Ga_m=Ga;


%% solution
for tt=2:length(Iem)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %solve the depth of upper ocean (DD) changes as the verturning changes
    %if overturning does not changes then DD is constant
    if tt==2
         DD(tt)=DD(tt-1)+Deltat*(VTek(tt-1)+VTed(tt-1)+VTna(tt-1)+VTu(tt-1))./Ar_gna;%Euler
    end
    %Adams-Bashforth 2nd order method for second time step
    if tt==3
        DD(tt)=DD(tt-1)+(1/Ar_gna)*Deltat*((3/2)*(VTek(tt-1)+VTed(tt-1)+VTna(tt-1)+VTu(tt-1))-(1/2)*(VTek(tt-2)+VTed(tt-2)+VTna(tt-2)+VTu(tt-2)));
    end
    %Adams-Bashforth 3rd order method for the rest time steps
    if tt>=4
         DD(tt)=DD(tt-1)+(1/Ar_gna)*Deltat*((23/12)*(VTek(tt-1)+VTed(tt-1)+VTna(tt-1)+VTu(tt-1))...
            -(4/3)*(VTek(tt-2)+VTed(tt-2)+VTna(tt-2)+VTu(tt-2))+(5/12)*(VTek(tt-3)+VTed(tt-3)+VTna(tt-3)+VTu(tt-3)));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update the terrestrial boxes
    Iveg(tt) = Iveg(tt-1) + dI_v(tt-1)*Deltat./(60^2*24*365);
    Isoil(tt) = Isoil(tt-1) + dI_s(tt-1)*Deltat./(60^2*24*365);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update the atmospheric CO2
    chi(tt)=chi(tt-1)-Deltat.*(Ar_o.*(Fa_all(tt-1)+Fa_all_pre))./Na...
        +(Iem(tt)-Iem(tt-1) - Iweath(tt) + Iweath(tt-1)).*10^15./(12.01*Na)...
        -(Iveg(tt)+Isoil(tt)-Iveg(tt-1)-Isoil(tt-1)).*10^15./(12.01*Na);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update radiative response and heat budget
    DR(tt)=alpha*(log(chi(tt)/chi(1))); 
    NTOA(tt)=(DR(tt)-lambda.*DTa(tt-1));
    Natm(tt)=CC.*(DTa(tt-1)-(T_cont(tt-1)-T_cont(1)));
    N(tt)=NTOA(tt)+Natm(tt);%N that goes into the ocean (almost all NTOA)
    
    %temperature evolution at each box (essentially solving generic tracer
    %equation)
    T_s(tt)=T_s(tt-1)+(Deltat./(Ar_s*D_s.*rho_o*cp_o)).*Hs(tt-1)+(Deltat./(Ar_s*D_s))*(VTRO(tt-1)*T_d(tt-1)-VTRO(tt-1)*T_s(tt-1));
    T_u(tt)=T_u(tt-1)*((DD(tt-1)-h_m)./(DD(tt)-h_m))+(Deltat./(Ar_gna*(DD(tt)-h_m)))*(FRAC*VTRO(tt-1)*T_s(tt-1)+(1-FRAC)*VTRO(tt-1)*T_m(tt-1)+VTna(tt-1)*T_u(tt-1)+VTu(tt-1)*T_d(tt-1));
    T_n(tt)=T_n(tt-1)+(Deltat./(Ar_n.*D_n.*rho_o*cp_o))*Hn_pre+(Deltat./(Ar_n*D_n.*rho_o*cp_o))*Hn(tt-1)+(Deltat./(Ar_n*D_n))*(-VTna(tt-1)*T_u(tt-1) +VTna(tt-1)*T_n(tt-1));
    T_m(tt)=T_m(tt-1)+(Deltat./(Ar_gna.*h_m.*rho_o*cp_o)).*Hm_pre+Deltat./(Ar_gna*h_m*rho_o*cp_o)*Hm(tt-1)+Deltat./(Ar_gna*h_m)*((1-FRAC)*VTRO(tt-1)*T_s(tt-1)-(1-FRAC)*VTRO(tt-1)*T_m(tt-1));
    T_d(tt)=T_d(tt-1)*((Ar_gna*(D-DD(tt-1))+Ar_s*(D-D_s)+Ar_n*(D-D_s))./(Ar_gna*(D-DD(tt))+Ar_s*(D-D_s)+Ar_n*(D-D_s)))...
              +(Deltat./(Ar_gna*(D-DD(tt))+Ar_s*(D-D_s)+Ar_n*(D-D_s)))*(-VTna(tt-1)*T_n(tt-1)-VTRO(tt-1)*T_d(tt-1)-VTu(tt-1)*T_d(tt-1));   
    
    T_cont(tt)=(T_s(tt).*Ar_s+T_n(tt).*Ar_n+T_m(tt).*Ar_gna)./Ar_o;
    % distribute where the anthopogenic heat goes
    % here I assume that the anthopogenic heat flux is not distributed
    % equally with 4/10 into the Souther ocean 2/10 into
    % tropics and subtropics and 4/10 into the north high latitudes
    % you can change as you choose
    Hm(tt)=Ar_o*N(tt)*0.2;
    Hn(tt)=Ar_o*N(tt)*0.4;
    Hs(tt)=Ar_o*N(tt)*0.4;
    
    DTa(tt)=(CC*(T_cont(tt)-T_cont(1))+rho_a*cp_a*h_a*DTa(tt-1)/Deltat)/(rho_a*cp_a*h_a/Deltat+CC);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update carbon flux and carbon budget
    
    %DIC evolution at each box (essentially solving generic tracer
    %equation)
    %16 July 2019- these have been updated with alkalinity inputs from HCO3
    %for the SH, NH, and ML boxes
    dic_s(tt)=dic_s(tt-1)+(Deltat./(Ar_s*D_s.*rho_o)).*(Fa_s(tt-1)+Ga_s(tt-1)).*Ar_s+(Deltat./(Ar_s*D_s))*(VTRO(tt-1)*dic_d(tt-1)-VTRO(tt-1)*dic_s(tt-1));
    dic_u(tt)=dic_u(tt-1)*((DD(tt-1)-h_m)./(DD(tt)-h_m))+...
        +(Deltat./(Ar_gna*(DD(tt)-h_m)))*(FRAC*VTRO(tt-1)*dic_s(tt-1)+(1-FRAC)*VTRO(tt-1)*dic_m(tt-1)+VTna(tt-1)*dic_u(tt-1)+VTu(tt-1)*dic_d(tt-1));
    dic_n(tt)=dic_n(tt-1)+(Deltat./(Ar_n*D_n.*rho_o))*(Fa_n(tt-1)+Fa_n_pre(end)+Ga_n(tt-1))*Ar_n+(Deltat./(Ar_n*D_n))*(-VTna(tt-1)*dic_u(tt-1)+...
         +VTna(tt-1)*dic_n(tt-1));
    dic_m(tt)=dic_m(tt-1)+Deltat./(Ar_gna*h_m*rho_o)*(Fa_m(tt-1)+Fa_m_pre(end)+Ga_m(tt-1))*Ar_gna+Deltat./(Ar_gna*h_m)*((1-FRAC)*VTRO(tt-1)*dic_s(tt-1)-(1-FRAC)*VTRO(tt-1)*dic_m(tt-1));
    dic_d(tt)=dic_d(tt-1)*((Ar_gna*(D-DD(tt-1))+Ar_s*(D-D_s)+Ar_n*(D-D_s))./(Ar_gna*(D-DD(tt))+Ar_s*(D-D_s)+Ar_n*(D-D_s)))...
              +(Deltat./(Ar_gna*(D-DD(tt))+Ar_s*(D-D_s)+Ar_n*(D-D_s)))*(-VTna(tt-1)*dic_n(tt-1)-VTRO(tt-1)*dic_d(tt-1)-VTu(tt-1)*dic_d(tt-1));
          
    %Alkalinity evolution at each box (essentially solving generic tracer
    %equation)
    %included 16 July 2019 --- these are the same equations as in the lines
    %above for carbon, but Ga's are converted to moleq/s
    %THIS IS SAME AS FOR GNA_F_1
    alk_s(tt)=alk_s(tt-1)+(Deltat./(Ar_s*D_s.*rho_o)).*(Ga_s(tt-1)).*Ar_s+(Deltat./(Ar_s*D_s))*(VTRO(tt-1)*alk_d(tt-1)-VTRO(tt-1)*alk_s(tt-1));
    alk_u(tt)=alk_u(tt-1)*((DD(tt-1)-h_m)./(DD(tt)-h_m))+...
        +(Deltat./(Ar_gna*(DD(tt)-h_m)))*(FRAC*VTRO(tt-1)*alk_s(tt-1)+(1-FRAC)*VTRO(tt-1)*alk_m(tt-1)+VTna(tt-1)*alk_u(tt-1)+VTu(tt-1)*alk_d(tt-1));
    alk_n(tt)=alk_n(tt-1)+(Deltat./(Ar_n*D_n.*rho_o))*(Ga_n(tt-1))*Ar_n+(Deltat./(Ar_n*D_n))*(-VTna(tt-1)*alk_u(tt-1)+...
         +VTna(tt-1)*alk_n(tt-1));
    alk_m(tt)=alk_m(tt-1)+Deltat./(Ar_gna*h_m*rho_o)*(Ga_m(tt-1))*Ar_gna+Deltat./(Ar_gna*h_m)*((1-FRAC)*VTRO(tt-1)*alk_s(tt-1)-(1-FRAC)*VTRO(tt-1)*alk_m(tt-1));
    alk_d(tt)=alk_d(tt-1)*((Ar_gna*(D-DD(tt-1))+Ar_s*(D-D_s)+Ar_n*(D-D_s))./(Ar_gna*(D-DD(tt))+Ar_s*(D-D_s)+Ar_n*(D-D_s)))...
              +(Deltat./(Ar_gna*(D-DD(tt))+Ar_s*(D-D_s)+Ar_n*(D-D_s)))*(-VTna(tt-1)*alk_n(tt-1)-VTRO(tt-1)*alk_d(tt-1)-VTu(tt-1)*alk_d(tt-1));
    
    % Estimate the anthopogenic ocean carbon exchange      
    dic_carb(tt)=(dic_s(tt).*Ar_s+dic_n(tt).*Ar_n+dic_m(tt).*Ar_gna)./Ar_o;
    At(tt)=(alk_s(tt).*Ar_s+alk_n(tt).*Ar_n+alk_m(tt).*Ar_gna)./Ar_o;
    [co2new H1 pH1,K1,K2,Ko]= co2_follows(T_cont(tt),Sref,dic_carb(tt),Pho,Sil,At(tt),pH_carb(tt-1));
    H_carb(tt)=H1;
    pH_carb(tt)=pH1;
    CO2_carb(tt)=dic_carb(tt)./(1+(K1/H_carb(tt))+(K1.*K2./(H_carb(tt)^2)));
    Fa_all(tt)=-rho_o.*Kg.*(CO2_carb(tt)-Ko.*P.*chi(tt))-Fa_all_pre(end); 
    % here I assume that the anthopogenic carbon flux is not distributed
    % equally with 4/10 into the Souther ocean 2/10 into
    % tropics and subtropics and 4/10 into the north high latitudes
    % you can change as you choose
    Fa_s(tt)=0.4*Ar_o*Fa_all(tt)/Ar_s;Fa_n(tt)=0.4*Ar_o*Fa_all(tt)/Ar_n;Fa_m(tt)=0.2*Ar_o*Fa_all(tt)/Ar_gna;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update the terrestrial system carbon fluxes using Phil's WASP code
    % citation: Goodwin 

    f1(tt) = 1 + dNPP_dT ./ NPP0 .* DTa(tt); %linear response of NPP to temperatures
    f2(tt) = 1 + gamma_co2 .* (log(chi(tt)) - log(chi_pre)); %empirical CO2 fertilisation
    f3(tt) = (Isoil0./NPP0) ./ (Isoil0./NPP0 + dtau_dT.*DTa(tt));

    dI_v(tt) = NPP0.*f1(tt).*f2(tt) - LL0.*Iveg(tt)./Iveg0;
    dI_s(tt) = LL0.*Iveg(tt)./Iveg0 - SR0.*f3(tt).*Isoil(tt)./Isoil0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update the volume transports (i.e., changes in circulation)
    %here the sinking in the north atlantic is a function of the heat uptake
    % and the warming gradients between upper and deep ocean at
    % tropics-subtropics (you can set other function if you want)
    VTna(tt)=VTna(1)+(Ar_gna*(N(tt)))/(rho_o*cp_o*(((DD(tt)-h_m)*T_u(tt)+h_m*T_m(tt))./DD(tt)-T_d(tt)));
    %IF YOU WANT TO KEEP CIRCULATION CONSTANT USE :VTna(tt)=VTna(1); 
    % the rest of the volume transports are updated based on changes in
    % thickness of upper ocean
    VTek(tt)=VTek(1);%the Ekman part does not change only the eddy
    VTu(tt)=VTuh/DD(tt);
    VTed(tt)=-VTedh*DD(tt);
    VTRO(tt)=VTek(tt)+VTed(tt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extra for analysis visualisation and validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate carbon inventories at each box in PgC
Ia=Na*chi*12.01/10^15;
Iu=dic_u.*((DD-h_m).*Ar_gna).*rho_o*12.01*10^(-15);
Id=(Ar_s*(D-D_s)+Ar_n*(D-D_n)+Ar_gna*(D-DD)).*dic_d*rho_o*12.01*10^(-15);
Is=dic_s.*(D_s)*Ar_s.*rho_o*12.01*10^(-15);
In=dic_n.*(D_n)*Ar_n.*rho_o*12.01*10^(-15);
It=dic_m.*(h_m*Ar_gna).*rho_o*12.01*10^(-15);

%check the budget closes Iocean+Ia=Iem
Iocean=It+Iu+Is+In+Id;
Iland=Iveg+Isoil;

%estimate heat content of the ocean it should matches the ocean heat uptake
Qupper=Ar_gna*rho_o*cp_o.*T_u.*(DD-h_m);
Qdeep=rho_o*cp_o.*T_d.*(Ar_gna.*(D-DD)+Ar_s.*(D-D_s)+Ar_n.*(D-D_n));
Qsouth=rho_o*cp_o.*T_s.*Ar_s.*(D_s);
Qnorth=rho_o*cp_o.*T_n.*Ar_n.*(D_n);
Qml=Ar_gna*rho_o*cp_o.*h_m.*(T_m);
Qtotalr=Qupper+Qdeep+Qsouth+Qnorth+Qml;
N_extra=diff(Qtotalr)./(Ar_o*Deltat);

MM=(1:length(Iem))/(365);%in years

clearvars -except DTa time lambda chi pH_carb tau FRAC Iem Rate_weath halftime Iocean Iland Ia