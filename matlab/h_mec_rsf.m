% ========================================
% H-MEC: Hydro-Mechanical Earthquake Cycle
% Computational Earthquake Physics
% ETH Zurich, 2022
% 
% Dal Zilio, L., Hegyi, B., Behr, W. M., & Gerya, T. (2022)
% Hydro-mechanical earthquake cycles in a 
% poro-visco-elasto-plastic fluid-bearing fault.
% DOI: http://arxiv.org/abs/2201.11786
% ========================================
% Solving of compressible Stokes+continuity equations
% with variable viscosity and Darsi+continuity equations
% Total X-Stokes: dSIGMAxxt'/dx+dSIGMAxyt'/dy-dPt/dx=-RHOt*gx
% Total Y-Stokes: dSIGMAyxt'/dx+dSIGMAyyt'/dy-dPt/dy=-RHOt*gy
% Solid Continuity: dVxs/dx+dVys/dy+(Pt-Pf)/ETAbulk=0
% Fluid X-Darsi: -ETAfluid/K*VxD-dPf/dx=-RHOf*gx
% Fluid Y-Darsi: -ETAfluid/K*VyD-dPf/dy=-RHOf*gy
% Fluid Continuity: dVxD/dx+dVyD/dy-(Pt-Pf)/ETAbulk=0
% + staggered grid
% + P-v formulation
% + advecting material fileds by markers
% ========================================
warning off
% ========================================
% Load mat file
fdata=fopen('file.txt','rt');
timestep=fscanf(fdata,'%d',1);
fclose(fdata);
if(timestep>0)
namemat    =  ['h_mec_',num2str(timestep)];
load(namemat);
timestep=timestep+1;
else
% ------------------------------
% ------------------------------
% 2) Define Numerical model
% Eulerian basic grid
xsize=5000.0; % size in horizontal direction, m
ysize=2000.0; % size in vertical direction, m
Nx=331;  % number of grid steps in horizontal directions
Ny=31; % number of grid steps in vertical direction
% Eulerian Staggered Grid
Nx1=Nx+1; % Number of horizontal lines for staggered grid
Ny1=Ny+1; % Number of vertical lines for staggered grid
% ========= 
seismic_cycles_data = 1;
timesum_plus        = 0;
line_fault          = (Ny-1)/2+1;
% =========
% Coordinates
dx=xsize/(Nx-1); % grid step in horizontal direction, m
dy=ysize/(Ny-1); % grid step in vertical direction, m
xbeg=0;
xend=xsize;
ybeg=0;
yend=ysize;
x=xbeg:dx:xend; % horizontal coordinates of basic grid points
y=ybeg:dy:yend; % vertical coordinates of basic grid points
xvx=xbeg:dx:xend+dx; % Horizontal coordinates of Vx-nodes
yvx=ybeg-dy/2:dy:yend+dy/2; % Vertical coordinates of Vx-nodes
xvy=xbeg-dx/2:dx:xend+dx/2; % Horizontal coordinates of Vy-nodes
yvy=ybeg:dy:yend+dy; % Vertical coordinates of Vy-nodes
xp=xbeg-dx/2:dx:xend+dx/2; % Horizontal coordinates of P-nodes
yp=ybeg-dy/2:dy:yend+dy/2; % Vertical coordinates of Vx-nodes

%          Bulk     Fault   
arsfm =   [0.020    0.007  ]; % a-parameter of RSF
brsfm =   [0.001    0.018  ]; % b-parameter of RSF
lrsfm =   [0.020    0.0010 ]; % L-parameter of RSF (characteristic slip distance)
omm   =   [15      -5      ]; % State

V0    =   1e-9;                   % Reference slip velocity of RSF, m/s
alpha =   29;                     % Scaling factor for shear viscosity of porous matrix,

% Brittle/plastic rheology
cohes=0.0e+6; %Cohesion, Pa
friction=0.6; % Internal friction coefficient confined
dilatation=0.00; % Dilatation coefficient confined
tensile=1.0; % Internal friction coefficient tensile
shearmod=3.0e+10; % Shear modulus
BETTAFLUID=1e-8;  % 4.0e-10; % Compressibility of fluid, 1/Pa
BETTASOLID=2e-11; % 2.5e-11; % Compressibility of solid, 1/Pa
PORDIL=0.15; % Upper limit of porosity for dilatation
ESLIP0=1e+04; % ETHAslip LN model, Pa s
VSLIP0=0.5e-9; % Characteristic velocity, m/s
faultwidth=dx; %250; % Characteristic fault width, m

% Constants
gx=0; % Horizontal gravity, m/s^2
gy=0; % Vertical gravity, m/s^2
SIGMALOAD=0; % SIGMAyy at the top, Pa
ETATOP=1e-3; % Viscosity of the sticky water, Pa*s

% Limits
pormin=1e-4; % Min porosity limit
pormax=1-pormin; % Max porosity limit
etamin=1e-3; % Lower shear viscosity cutoff
etamax=1e+50; % Upper shear viscosity cutoff
etabmin=1e-3; % Lower bulk viscosity cutoff
etabmax=1e+25; % Upper bulk viscosity cutoff
gmin=1e+08; % Lower shear modulus cutoff
gmax=1e+11; % Upper shear modulus cutoff
kkkmin=1e-22; % Lower Darsi viscosity cutoff
kkkmax=1e-12; % Upper Darsi viscosity cutoff
kdtmin=0.99; % Max change in porosity per timestep


% Limits
dsubgrids=0; % Subgrid diffusion for stresses
dsubgridv=0; % Subgrid diffusion for velocity
stpmax=1e-1;%/dy*faultwidth; % Max gridstep fraction for marker displacement in the channel
stpmax1=1e-2;%/dy*faultwidth; % Max gridstep fraction for marker displacement
dvmax=0.5; % Max relative velocity change per timestep


OM0=omm(1)*ones(Ny,Nx); % Old state parameter
OM=omm(1)*ones(Ny,Nx); % State parameter
ARSF=arsfm(1)*ones(Ny,Nx); % a-parameter of RSF
BRSF=brsfm(1)*ones(Ny,Nx); % b-parameter of RSF
LRSF=lrsfm(1)*ones(Ny,Nx); % L-parameter of RSF

% Unknown parameters
pt=zeros(Ny1,Nx1); % Total pressure
vxs=zeros(Ny1,Nx1); % Solid vx-velocity
vys=zeros(Ny1,Nx1); % Solid vy-velocity
vzs=zeros(Ny1,Nx1); % Solid vy-velocity
pf=zeros(Ny1,Nx1); % Fluid pressure
vxD=zeros(Ny1,Nx1); % Darsi vx-velocity
vyD=zeros(Ny1,Nx1); % Darsi vy-velocity

% Nodal arrays
% Basic nodes
RHO=zeros(Ny,Nx);
ETA=zeros(Ny,Nx);
ETA0=zeros(Ny,Nx);
IETAPLB=zeros(Ny,Nx);
SXY=zeros(Ny,Nx);
SXY0=zeros(Ny,Nx);

% FABIAN
SZX0=zeros(Ny,Nx);
SZY0=zeros(Ny,Nx);
% ==============
YNY0=zeros(Ny,Nx);
KKK=zeros(Ny,Nx);
GGG=zeros(Ny,Nx);
COHC=zeros(Ny,Nx);
COHT=zeros(Ny,Nx);
FRIC=zeros(Ny,Nx);
FRIT=zeros(Ny,Nx);
DILC=zeros(Ny,Nx);
TTT=zeros(Ny,Nx);
EIIB=zeros(Ny,Nx);
STRPLB=zeros(Ny,Nx);
VSLIPB=zeros(Ny,Nx);
% ===========================
% Frictional parameters

bulk_up=line_fault-1; % right side of the bulk
bulk_lw=line_fault+1; % left side of the bulk

a_RSF_VS_UP=arsfm(1);
a_RSF_VS_BT=arsfm(1);
a_RSF_VW=arsfm(2);

b_RSF_VS=brsfm(1);
b_RSF_VW=brsfm(2);

BRSF(:,:)=b_RSF_VW;
ARSF(:,:)=a_RSF_VW;

L_VS = lrsfm(1);
L_VW = lrsfm(2);

LRSF(:,:)=L_VS;
% ===========================
% transition between VW to VS

TS_1 = 100;
TS_2 = 450;

TS_3 = xsize-450;
TS_4 = xsize-100;

% Define Fault
for j=1:1:Nx
    for i=1:1:Ny

        OM0(i,j)=omm(1);

        if (i==line_fault) 
            OM0(i,j)=omm(2);
        end

        if(x(j)<TS_1 && i==line_fault) 
            BRSF(i,j)=b_RSF_VS;
            ARSF(i,j)=a_RSF_VS_UP;
            LRSF(i,j)=L_VS;
        end

        if(x(j)>=TS_1 && x(j)<TS_2 && i==line_fault)
            BRSF(i,j)=b_RSF_VS-(b_RSF_VS-b_RSF_VW)*((x(j)-TS_1)/(TS_2-TS_1));
            ARSF(i,j)=a_RSF_VS_UP-(a_RSF_VS_UP-a_RSF_VW)*((x(j)-TS_1)/(TS_2-TS_1));
            LRSF(i,j)=L_VS-(L_VS-L_VW)*((x(j)-TS_1)/(TS_2-TS_1));
        end

        if(x(j)>=TS_2 && x(j)<=TS_3 && i==line_fault)
            BRSF(i,j)=b_RSF_VW;
            ARSF(i,j)=a_RSF_VW;
            LRSF(i,j)=L_VW;
        end

        if(x(j)>TS_3 && x(j)<=TS_4 && i==line_fault) 
            BRSF(i,j)=b_RSF_VW-(b_RSF_VW-b_RSF_VS)*((x(j)-TS_3)/(TS_4-TS_3));
            ARSF(i,j)=a_RSF_VW-(a_RSF_VW-a_RSF_VS_BT)*((x(j)-TS_3)/(TS_4-TS_3));
            LRSF(i,j)=L_VW-(L_VW-L_VS)*((x(j)-TS_3)/(TS_4-TS_3));
        end
        if(x(j)>TS_4 && i==line_fault)
            BRSF(i,j)=b_RSF_VS;
            ARSF(i,j)=a_RSF_VS_BT;
            LRSF(i,j)=L_VS;
        end
    end
end

x_vs_beg = TS_2+200;
checkpt=0;

KKK_barrier=1e-22;
KKK(:,:)=1e-20;

for j=1:1:Nx
    for i=1:1:Ny    
        
        if x(j)>=x_vs_beg && x(j)>TS_2 && x(j)<TS_3
            
            % Min and Max width barriers
            if checkpt==0
                x_vs_end = (200-90).*rand(1,1)+90;
                checkpt=1;
            end
            
            if (x(j)>=x_vs_beg && x(j)<=(x_vs_beg+x_vs_end))
                KKK(i,j)=KKK_barrier;
            end
            
            % Min and Max width seismic patch --> compare it to h* 
            if x(j)>(x_vs_beg+x_vs_end)
                x_vs_beg =x_vs_beg+x_vs_end+(250-180).*rand(1,1)+180;
                checkpt=0;
            end
        end
    end
end

p_diff          = 30e6;
p_fluid         = 10e6;
p_fluid_barrier = 35e6;
PCONF           = zeros(Ny1,Nx1)+p_fluid;
PTFDIFF         = zeros(Ny1,Nx1)+p_diff;

for i=1:1:Ny
    for j=1:1:Nx
        if(KKK(i,j)<=KKK_barrier)% || x(j)<=TS_1 || x(j)>=TS_4)
            PCONF(i,j)         = p_fluid_barrier;
            PTFDIFF(i,j)       = p_diff-(p_fluid_barrier-p_fluid);
	    ARSF(line_fault,j) = a_RSF_VS_UP; 
        else
            PCONF(i,j)=p_fluid;
            PTFDIFF(i,j)=p_diff;
        end
    end
end

ARSF(line_fault,:) = smoothdata(ARSF(line_fault,:),'gaussian',4);
KKK_smooth         = smoothdata(KKK(line_fault,:),'gaussian',4);
Pt_smooth          = smoothdata(PTFDIFF(line_fault,:),'gaussian',4);
Pf_smooth          = smoothdata(PCONF(line_fault,:),'gaussian',4);

for i=1:1:Ny
    PTFDIFF(i,:) = Pt_smooth;
    PCONF(i,:)   = Pf_smooth;
    KKK(i,:)     = KKK_smooth;
end

PT0=(PCONF+PTFDIFF);
PF0=PCONF;

por       = 0.01; % porosity
POR0      = 0.01; % reference porosity
POR       = zeros(Ny1,Nx1)+por;
KKK       = KKK*(por/POR0)^3;
eta_fluid = 1e-3; % Fluid viscosity
ETAD      = eta_fluid./KKK; % Darsi "viscosity"
% =================
% nucleation size
hstar=min(pi/2*shearmod*lrsfm(2).*brsfm(2)./(brsfm(2)-arsfm(2)).^2./max(max(PTFDIFF)));
% cohesive zone size
coh = min(9/32*pi*shearmod*lrsfm(2)./brsfm(2)./max(max(PTFDIFF)));

% Print information about discretization
fprintf('>> Grid size = %.2f (m)\n', dx);
fprintf('>> VW width = %.2f (km)\n', (TS_4+TS_3)/2/1e3);
fprintf('>> Critical nucleation size = %.2f (m)\n',hstar);
fprintf('>> Cohesive zone = %.2f (m)\n',coh);
fprintf('>> Cohesive zone/dx = %.2f \n',coh/dx);

OM=OM0;

% Pressure nodes
ETAB=zeros(Ny1,Nx1);
ETAB0=zeros(Ny1,Nx1);
ETAP=zeros(Ny1,Nx1);
ETAP0=zeros(Ny1,Nx1);
GGGP=zeros(Ny1,Nx1);
GGGB=zeros(Ny1,Nx1);
PTF0=zeros(Ny1,Nx1);
pt_ave=zeros(Ny1,Nx1);
pf_ave=zeros(Ny1,Nx1);
SXX=zeros(Ny1,Nx1);
SXX0=zeros(Ny1,Nx1);
SYY=zeros(Ny1,Nx1);
SYY0=zeros(Ny1,Nx1);
DILP=zeros(Ny1,Nx1);
% Vx nodes
RHOX=zeros(Ny1,Nx1);
RHOFX=zeros(Ny1,Nx1);
VX0=zeros(Ny1,Nx1);
VXF0=zeros(Ny1,Nx1);
% Vy nodes
RHOY=zeros(Ny1,Nx1);
RHOFY=zeros(Ny1,Nx1);
VY0=zeros(Ny1,Nx1);
VYF0=zeros(Ny1,Nx1);

% Vz nodes
VZ0=zeros(Ny1,Nx1);

% Lagrangian solid markers
Nxm=(Nx-1)*3; % Marker resolution in x-dection
Nym=(Ny-1)*4; % Marker resolution in y direction
dxms=xsize/Nxm; % Standard marker horizontal step
dyms=ysize/Nym; % Standard marker vertical step
marknum=Nxm*Nym; % Total number of markers
rhom=zeros(marknum,1); % Density of solid
etasm=zeros(marknum,1); % Standard shear viscosity of bulk
etam=zeros(marknum,1); % Shear viscosity of bulk
etabm=zeros(marknum,1); % Bulk viscosity of bulk
cohescm=zeros(marknum,1); % Cohesion for confined fracture of solid
frictcm=zeros(marknum,1); % friction for confined fracture of solid
dilatcm=zeros(marknum,1); % dilatation for confined fracture of solid
cohestm=zeros(marknum,1); % Cohesion for tensile fracture of solid
fricttm=zeros(marknum,1); % friction for tensile fracture of solid
rhofm=zeros(marknum,1); % Density of fluid
tm=zeros(marknum,1); % Marker rock type
xm=zeros(marknum,1); % Horizontal coordinates of solid markers
ym=zeros(marknum,1); % Vertical coordinates of solid markers
sxxm=zeros(marknum,1); % Marker SIGMAxx', Pa
syym=zeros(marknum,1); % Marker SIGMAyy', Pa
sxym=zeros(marknum,1); % Marker SIGMAxy', Pa
szxm=zeros(marknum,1); % Marker SIGMAzx', Pa
szym=zeros(marknum,1); % Marker SIGMAzy', Pa
gsm=zeros(marknum,1); % Standard shear modulus of bulk, Pa
gm=zeros(marknum,1); % Shear modulus of bulk, Pa
vx0m=zeros(marknum,1); % Marker horizontal velocity
vy0m=zeros(marknum,1); % Marker vertical velocity
vz0m=zeros(marknum,1); % Marker vertical velocity
ptfm=zeros(marknum,1); % Pt-Pf, Pa
amursfm=zeros(marknum,1); % RSF a/mu parameter

m=1;
for jm=1:1:Nxm
    for im=1:1:Nym
        % Define randomized regular coordinates
        xm(m)=xbeg+dxms/2+(jm-1)*dxms+(rand-0.5)*dxms;
        ym(m)=ybeg+dyms/2+(im-1)*dyms+(rand-0.5)*dyms;
        %Matrix
        tm(m)=1;
        rhom(m)=3000;
        etasm(m)=10^21;
        gsm(m)=shearmod;
        cohescm(m)=cohes;
        cohestm(m)=cohes;
        frictcm(m)=friction;
        dilatcm(m)=dilatation;
        fricttm(m)=tensile;
        rhofm(m)=1000;
        etam(m)=etasm(m)*exp(-alpha*por);
        etabm(m)=etam(m)/por;%/(1-por);
        gm(m)=gsm(m);%*(1-por);

        % Air, wedge, slab
        if(ym(m)<y(bulk_up) || ym(m)>y(bulk_lw))
            tm(m)=-1;
            etam(m)=1e+25;
            etasm(m)=1e+25;
            rhom(m)=3000;
            cohescm(m)=cohes*1e+3;
            cohestm(m)=cohes*1e+3;
            frictcm(m)=1.0;%friction;
            dilatcm(m)=dilatation;
            fricttm(m)=tensile;
        end

        % Update marker index
        m=m+1;
    end
end

% 3) Defining global matrixes
% according to the global number of unknowns
N=Nx1*Ny1*6; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients in the left part
% L=spalloc(N,N,N*10); % Matrix of coefficients in the left part
R=zeros(N,1); % Vector of the right parts of equations


% 4) Boundary conditions
bcleft=1;
bcright=1;
bcupper=-1e-9;
bclower=1e-9;
bcvyflower=0;%-1e-12;
% bcvxfleft=bcvyflower*xsize/ysize;


% Entering time cycle
timesum=0;
dtelastic0=3e+7; % elastic timestep
dt=1e4;
dtmin=1e-4;

ascale=1e+0;
dtreset=1e-20; % Minimum timestep for velocity reset

% Timesteps between visualization frames
savematstep=300;  % .mat storage periodicity
nname ='h_mec_';  % mat filename
niterglobal=10000; % Max number of global iterations
ynlastmax=499; % Max number of iterations without changing timestep
dtstep=5; % Min number of plastic iterations before changing dt
derrdstp=-0.000; % Min error decrease per step
vratiomax=0.001; % Max change of velocity between iterations
kdtetmin=0.01; % Herrendoerfer 2017 timestepping constant
dtkoef=1.1; % Koefficient for decrease of previous timestep
dtkoefup=1.2; % Koefficient for increase of previous timestep
dtkoefv=1.001; % koef for velocity-based timestep reduction
errmin=1e+3; % min LSQ err for stoping iterations
etawt=0.0;
syieldmin=1e3;
tyield=1;
timestep=1;
maxvxsmod=0;
minvxsmod=0;
maxvysmod=0;
minvysmod=0;
maxvzsmod=0;
minvzsmod=0;
yndtdecrease=1;

end


for timestep=timestep:1:1000000
    
    
    % Interpolate ETA, RHO to nodal points
    % Basic nodes
    RHOSUM=zeros(Ny,Nx);
    ETASUM=zeros(Ny,Nx);
    TTTSUM=zeros(Ny,Nx);
    SXYSUM=zeros(Ny,Nx);
    GGGSUM=zeros(Ny,Nx);
    ETA0SUM=zeros(Ny,Nx);
    COHCSUM=zeros(Ny,Nx);
    FRICSUM=zeros(Ny,Nx);
    DILCSUM=zeros(Ny,Nx);
    COHTSUM=zeros(Ny,Nx);
    FRITSUM=zeros(Ny,Nx);
    WTSUM=zeros(Ny,Nx);
    
    OM0SUM=zeros(Ny,Nx); % Old state parameter
    OMSUM=zeros(Ny,Nx); % State parameter
    ARSFSUM=zeros(Ny,Nx); % a-parameter of RSF
    BRSFSUM=zeros(Ny,Nx); % b-parameter of RSF
    LRSFSUM=zeros(Ny,Nx); % L-parameter of RSF
    
    % Pressure nodes
    ETAPSUM=zeros(Ny1,Nx1);
    ETAP0SUM=zeros(Ny1,Nx1);
    ETAB0SUM=zeros(Ny1,Nx1);
    SXXSUM=zeros(Ny1,Nx1);
    SYYSUM=zeros(Ny1,Nx1);
    GGGPSUM=zeros(Ny1,Nx1);
    WTPSUM=zeros(Ny1,Nx1);
    % Vx nodes
    RHOXSUM=zeros(Ny1,Nx1);
    RHOFXSUM=zeros(Ny1,Nx1);
    VX0SUM=zeros(Ny1,Nx1);
    WTXSUM=zeros(Ny1,Nx1);
    % Vy nodes
    RHOYSUM=zeros(Ny1,Nx1);
    RHOFYSUM=zeros(Ny1,Nx1);
    VY0SUM=zeros(Ny1,Nx1);
    WTYSUM=zeros(Ny1,Nx1);
    
    VZ0SUM=zeros(Ny1,Nx1);
    
    % Cycle on markers
    for m=1:1:marknum
        
        % Cohesion, friction of porous matrix
        
        % Viscosity of porous matrix
        if(tm(m)~=0)
            cohescmm=cohescm(m)*(1-por);%*exp(-alpha*por);
            cohestmm=cohestm(m)*(1-por);%*exp(-alpha*por);
            frictcmm=frictcm(m);
            dilatcmm=dilatcm(m);
            fricttmm=fricttm(m);
            etasmm0=etasm(m);
            etamm0=etasm(m)*exp(-alpha*por);
            etamm=etam(m);
            % total density
            rhomm=rhom(m)*(1-por)+rhofm(m)*por;
        else
            cohescmm=cohescm(m);
            cohestmm=cohestm(m);
            frictcmm=frictcm(m);
            dilatcmm=0;
            fricttmm=fricttm(m);
            etasmm0=etasm(m);
            etamm0=etasm(m);
            etamm=etam(m);
            % total density
            rhomm=rhom(m);
        end
        % Matrix viscosity
        if(etamm0<etamin)
            etamm0=etamin;
        end
        if(etamm0>etamax)
            etamm0=etamax;
        end
        % Effective viscosity
        if(etamm<etamin)
            etamm=etamin;
        end
        if(etamm>etamax)
            etamm=etamax;
        end
        
        
        % Interpolate to basic nodes
        % [i,j]--------[i,j+1]
        %   |                |
        %   |    o m         |
        %   |                |
        % [i+1,j]-------[i+1,j+1]
        % Indexes and distances
        j=fix(xm(m)/dx)+1;
        i=fix(ym(m)/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx-1)
            j=Nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>Ny-1)
            i=Ny-1;
        end
        dxm=(xm(m)-x(j))/dx;
        dym=(ym(m)-y(i))/dy;
        
        
        % Updating
        % i,j
        wtmij=(1-dxm)*(1-dym);
        ETASUM(i,j)=ETASUM(i,j)+etamm*wtmij;
        VZ0SUM(i,j)=VZ0SUM(i,j)+vz0m(m)*rhomm*wtmij;
        RHOSUM(i,j)=RHOSUM(i,j)+rhomm*wtmij;
        TTTSUM(i,j)=TTTSUM(i,j)+tm(m)*wtmij;
        SXYSUM(i,j)=SXYSUM(i,j)+sxym(m)*wtmij;
        GGGSUM(i,j)=GGGSUM(i,j)+1/gm(m)*wtmij;
        ETA0SUM(i,j)=ETA0SUM(i,j)+etamm0*wtmij;
        COHTSUM(i,j)=COHTSUM(i,j)+cohestmm*wtmij;
        FRITSUM(i,j)=FRITSUM(i,j)+fricttmm*wtmij;
        COHCSUM(i,j)=COHCSUM(i,j)+cohescmm*wtmij;
        FRICSUM(i,j)=FRICSUM(i,j)+frictcmm*wtmij;
        DILCSUM(i,j)=DILCSUM(i,j)+dilatcmm*wtmij;
        
        
        WTSUM(i,j)=WTSUM(i,j)+wtmij;
        
        % i+1,j
        wtmi1j=(1-dxm)*(dym);
        ETASUM(i+1,j)=ETASUM(i+1,j)+etamm*wtmi1j;
        VZ0SUM(i,j)=VZ0SUM(i,j)+vz0m(m)*rhomm*wtmi1j;
        RHOSUM(i+1,j)=RHOSUM(i+1,j)+rhomm*wtmi1j;
        TTTSUM(i+1,j)=TTTSUM(i+1,j)+tm(m)*wtmi1j;
        SXYSUM(i+1,j)=SXYSUM(i+1,j)+sxym(m)*wtmi1j;
        GGGSUM(i+1,j)=GGGSUM(i+1,j)+1/gm(m)*wtmi1j;
        ETA0SUM(i+1,j)=ETA0SUM(i+1,j)+etamm0*wtmi1j;
        COHTSUM(i+1,j)=COHTSUM(i+1,j)+cohestmm*wtmi1j;
        FRITSUM(i+1,j)=FRITSUM(i+1,j)+fricttmm*wtmi1j;
        COHCSUM(i+1,j)=COHCSUM(i+1,j)+cohescmm*wtmi1j;
        FRICSUM(i+1,j)=FRICSUM(i+1,j)+frictcmm*wtmi1j;
        DILCSUM(i+1,j)=DILCSUM(i+1,j)+dilatcmm*wtmi1j;
        
        WTSUM(i+1,j)=WTSUM(i+1,j)+wtmi1j;
        
        
        % i,j+1
        wtmij1=(dxm)*(1-dym);
        ETASUM(i,j+1)=ETASUM(i,j+1)+etamm*wtmij1;
        VZ0SUM(i,j)=VZ0SUM(i,j)+vz0m(m)*rhomm*wtmij1;
        RHOSUM(i,j+1)=RHOSUM(i,j+1)+rhomm*wtmij1;
        TTTSUM(i,j+1)=TTTSUM(i,j+1)+tm(m)*wtmij1;
        SXYSUM(i,j+1)=SXYSUM(i,j+1)+sxym(m)*wtmij1;
        GGGSUM(i,j+1)=GGGSUM(i,j+1)+1/gm(m)*wtmij1;
        ETA0SUM(i,j+1)=ETA0SUM(i,j+1)+etamm0*wtmij1;
        COHTSUM(i,j+1)=COHTSUM(i,j+1)+cohestmm*wtmij1;
        FRITSUM(i,j+1)=FRITSUM(i,j+1)+fricttmm*wtmij1;
        COHCSUM(i,j+1)=COHCSUM(i,j+1)+cohescmm*wtmij1;
        FRICSUM(i,j+1)=FRICSUM(i,j+1)+frictcmm*wtmij1;
        DILCSUM(i,j+1)=DILCSUM(i,j+1)+dilatcmm*wtmij1;
        
        WTSUM(i,j+1)=WTSUM(i,j+1)+wtmij1;
        
        
        % i+1,j+1
        wtmi1j1=(dxm)*(dym);
        ETASUM(i+1,j+1)=ETASUM(i+1,j+1)+etamm*wtmi1j1;
        VZ0SUM(i,j)=VZ0SUM(i,j)+vz0m(m)*rhomm*wtmi1j1;
        RHOSUM(i+1,j+1)=RHOSUM(i+1,j+1)+rhomm*wtmi1j1;
        TTTSUM(i+1,j+1)=TTTSUM(i+1,j+1)+tm(m)*wtmi1j1;
        SXYSUM(i+1,j+1)=SXYSUM(i+1,j+1)+sxym(m)*wtmi1j1;
        GGGSUM(i+1,j+1)=GGGSUM(i+1,j+1)+1/gm(m)*wtmi1j1;
        ETA0SUM(i+1,j+1)=ETA0SUM(i+1,j+1)+etamm0*wtmi1j1;
        COHTSUM(i+1,j+1)=COHTSUM(i+1,j+1)+cohestmm*wtmi1j1;
        FRITSUM(i+1,j+1)=FRITSUM(i+1,j+1)+fricttmm*wtmi1j1;
        COHCSUM(i+1,j+1)=COHCSUM(i+1,j+1)+cohescmm*wtmi1j1;
        FRICSUM(i+1,j+1)=FRICSUM(i+1,j+1)+frictcmm*wtmi1j1;
        DILCSUM(i+1,j+1)=DILCSUM(i+1,j+1)+dilatcmm*wtmi1j1;
        
        WTSUM(i+1,j+1)=WTSUM(i+1,j+1)+wtmi1j1;
        
        % Interpolate to pressure nodes
        % [i,j]--------[i,j+1]
        %   |                |
        %   |    o m         |
        %   |                |
        % [i+1,j]-------[i+1,j+1]
        % Indexes and distances
        j=fix((xm(m)+dx/2)/dx)+1;
        i=fix((ym(m)+dy/2)/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        dxm=(xm(m)-xp(j))/dx;
        dym=(ym(m)-yp(i))/dy;
        
        % Updating
        % i,j
        
        wtmij=(1-dxm)*(1-dym);
        ETAPSUM(i,j)=ETAPSUM(i,j)+etamm*wtmij;
        SXXSUM(i,j)=SXXSUM(i,j)+sxxm(m)*wtmij;
        SYYSUM(i,j)=SYYSUM(i,j)+syym(m)*wtmij;
        GGGPSUM(i,j)=GGGPSUM(i,j)+1/gm(m)*wtmij;
        ETAP0SUM(i,j)=ETAP0SUM(i,j)+etamm0*wtmij;
        ETAB0SUM(i,j)=ETAB0SUM(i,j)+etasmm0*wtmij;
        WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
        
        % i+1,j
        
        wtmi1j=(1-dxm)*(dym);
        ETAPSUM(i+1,j)=ETAPSUM(i+1,j)+etamm*wtmi1j;
        SXXSUM(i+1,j)=SXXSUM(i+1,j)+sxxm(m)*wtmi1j;
        SYYSUM(i+1,j)=SYYSUM(i+1,j)+syym(m)*wtmi1j;
        GGGPSUM(i+1,j)=GGGPSUM(i+1,j)+1/gm(m)*wtmi1j;
        ETAP0SUM(i+1,j)=ETAP0SUM(i+1,j)+etamm0*wtmi1j;
        ETAB0SUM(i+1,j)=ETAB0SUM(i+1,j)+etasmm0*wtmi1j;
        WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
        
        % i,j+1
        
        wtmij1=(dxm)*(1-dym);
        ETAPSUM(i,j+1)=ETAPSUM(i,j+1)+etamm*wtmij1;
        SXXSUM(i,j+1)=SXXSUM(i,j+1)+sxxm(m)*wtmij1;
        SYYSUM(i,j+1)=SYYSUM(i,j+1)+syym(m)*wtmij1;
        GGGPSUM(i,j+1)=GGGPSUM(i,j+1)+1/gm(m)*wtmij1;
        ETAP0SUM(i,j+1)=ETAP0SUM(i,j+1)+etamm0*wtmij1;
        ETAB0SUM(i,j+1)=ETAB0SUM(i,j+1)+etasmm0*wtmij1;
        WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
        
        % i+1,j+1
        
        wtmi1j1=(dxm)*(dym);
        ETAPSUM(i+1,j+1)=ETAPSUM(i+1,j+1)+etamm*wtmi1j1;
        SXXSUM(i+1,j+1)=SXXSUM(i+1,j+1)+sxxm(m)*wtmi1j1;
        SYYSUM(i+1,j+1)=SYYSUM(i+1,j+1)+syym(m)*wtmi1j1;
        GGGPSUM(i+1,j+1)=GGGPSUM(i+1,j+1)+1/gm(m)*wtmi1j1;
        ETAP0SUM(i+1,j+1)=ETAP0SUM(i+1,j+1)+etamm0*wtmi1j1;
        ETAB0SUM(i+1,j+1)=ETAB0SUM(i+1,j+1)+etasmm0*wtmi1j1;
        WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;
        
        % Interpolate to Vx nodes
        % [i,j]--------[i,j+1]
        %   |                |
        %   |    o m         |
        %   |                |
        % [i+1,j]-------[i+1,j+1]
        % Indexes and distances
        j=fix((xm(m))/dx)+1;
        i=fix((ym(m)+dy/2)/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx-1)
            j=Nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        dxm=(xm(m)-xvx(j))/dx;
        dym=(ym(m)-yvx(i))/dy;
        
        % Updating
        % i,j
        wtmij=(1-dxm)*(1-dym);
        RHOXSUM(i,j)=RHOXSUM(i,j)+rhomm*wtmij;
        RHOFXSUM(i,j)=RHOFXSUM(i,j)+rhofm(m)*wtmij;
        VX0SUM(i,j)=VX0SUM(i,j)+vx0m(m)*rhomm*wtmij;
        WTXSUM(i,j)=WTXSUM(i,j)+wtmij;
        
        % i+1,j
        wtmi1j=(1-dxm)*(dym);
        RHOXSUM(i+1,j)=RHOXSUM(i+1,j)+rhomm*wtmi1j;
        RHOFXSUM(i+1,j)=RHOFXSUM(i+1,j)+rhofm(m)*wtmi1j;
        VX0SUM(i+1,j)=VX0SUM(i+1,j)+vx0m(m)*rhomm*wtmi1j;
        WTXSUM(i+1,j)=WTXSUM(i+1,j)+wtmi1j;
        
        % i,j+1
        wtmij1=(dxm)*(1-dym);
        RHOXSUM(i,j+1)=RHOXSUM(i,j+1)+rhomm*wtmij1;
        RHOFXSUM(i,j+1)=RHOFXSUM(i,j+1)+rhofm(m)*wtmij1;
        VX0SUM(i,j+1)=VX0SUM(i,j+1)+vx0m(m)*rhomm*wtmij1;
        WTXSUM(i,j+1)=WTXSUM(i,j+1)+wtmij1;
        
        % i+1,j+1
        wtmi1j1=(dxm)*(dym);
        RHOXSUM(i+1,j+1)=RHOXSUM(i+1,j+1)+rhomm*wtmi1j1;
        RHOFXSUM(i+1,j+1)=RHOFXSUM(i+1,j+1)+rhofm(m)*wtmi1j1;
        VX0SUM(i+1,j+1)=VX0SUM(i+1,j+1)+vx0m(m)*rhomm*wtmi1j1;
        WTXSUM(i+1,j+1)=WTXSUM(i+1,j+1)+wtmi1j1;
        
        
        % Interpolate to Vy nodes
        % [i,j]--------[i,j+1]
        %   |                |
        %   |    o m         |
        %   |                |
        % [i+1,j]-------[i+1,j+1]
        % Indexes and distances
        j=fix((xm(m)+dx/2)/dx)+1;
        i=fix((ym(m))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny-1)
            i=Ny-1;
        end
        dxm=(xm(m)-xvy(j))/dx;
        dym=(ym(m)-yvy(i))/dy;
        % Updating
        % i,j
        %     if(dxm<=0.5 && dym<=0.5)
        wtmij=(1-dxm)*(1-dym);
        RHOYSUM(i,j)=RHOYSUM(i,j)+rhomm*wtmij;
        RHOFYSUM(i,j)=RHOFYSUM(i,j)+rhofm(m)*wtmij;
        VY0SUM(i,j)=VY0SUM(i,j)+vy0m(m)*rhomm*wtmij;
        WTYSUM(i,j)=WTYSUM(i,j)+wtmij;
        %     end
        % i+1,j
        %     if(dxm<=0.5 && dym>=0.5)
        wtmi1j=(1-dxm)*(dym);
        RHOYSUM(i+1,j)=RHOYSUM(i+1,j)+rhomm*wtmi1j;
        RHOFYSUM(i+1,j)=RHOFYSUM(i+1,j)+rhofm(m)*wtmi1j;
        VY0SUM(i+1,j)=VY0SUM(i+1,j)+vy0m(m)*rhomm*wtmi1j;
        WTYSUM(i+1,j)=WTYSUM(i+1,j)+wtmi1j;
        %     end
        % i,j+1
        %     if(dxm>=0.5 && dym<=0.5)
        wtmij1=(dxm)*(1-dym);
        RHOYSUM(i,j+1)=RHOYSUM(i,j+1)+rhomm*wtmij1;
        RHOFYSUM(i,j+1)=RHOFYSUM(i,j+1)+rhofm(m)*wtmij1;
        VY0SUM(i,j+1)=VY0SUM(i,j+1)+vy0m(m)*rhomm*wtmij1;
        WTYSUM(i,j+1)=WTYSUM(i,j+1)+wtmij1;
        %     end
        % i+1,j+1
        %     if(dxm>=0.5 && dym>=0.5)
        wtmi1j1=(dxm)*(dym);
        RHOYSUM(i+1,j+1)=RHOYSUM(i+1,j+1)+rhomm*wtmi1j1;
        RHOFYSUM(i+1,j+1)=RHOFYSUM(i+1,j+1)+rhofm(m)*wtmi1j1;
        VY0SUM(i+1,j+1)=VY0SUM(i+1,j+1)+vy0m(m)*rhomm*wtmi1j1;
        WTYSUM(i+1,j+1)=WTYSUM(i+1,j+1)+wtmi1j1;
        %     end
        
    end
    
    % Computing ETA and RHO
    % Basic nodes
    for i=1:1:Ny
        for j=1:1:Nx
            if(WTSUM(i,j)>0)
                ETA(i,j)=ETASUM(i,j)/WTSUM(i,j);
                VZ0(i,j)=VZ0SUM(i,j)/WTSUM(i,j);
                RHO(i,j)=RHOSUM(i,j)/WTSUM(i,j);
                TTT(i,j)=(TTTSUM(i,j)/WTSUM(i,j));
                GGG(i,j)=1/(GGGSUM(i,j)/WTSUM(i,j));
                ETA0(i,j)=ETA0SUM(i,j)/WTSUM(i,j);
                COHT(i,j)=(COHTSUM(i,j)/WTSUM(i,j));
                FRIT(i,j)=(FRITSUM(i,j)/WTSUM(i,j));
                COHC(i,j)=(COHCSUM(i,j)/WTSUM(i,j));
                FRIC(i,j)=(FRICSUM(i,j)/WTSUM(i,j));
                DILC(i,j)=(DILCSUM(i,j)/WTSUM(i,j));
            end
        end
    end
    
    
    %Pressure nodes
    for i=1:1:Ny1
        for j=1:1:Nx1
            if(WTPSUM(i,j)>0)
                ETAP(i,j)=ETAPSUM(i,j)/WTPSUM(i,j);
                ETAB0(i,j)=ETAB0SUM(i,j)/WTPSUM(i,j)/POR(i,j);
                ETAP0(i,j)=ETAP0SUM(i,j)/WTPSUM(i,j);
                GGGP(i,j)=1/(GGGPSUM(i,j)/WTPSUM(i,j));
            end
        end
    end
    % Vx nodes
    for i=1:1:Ny1
        for j=1:1:Nx
            if(WTXSUM(i,j)>0)
                RHOX(i,j)=RHOXSUM(i,j)/WTXSUM(i,j);
                RHOFX(i,j)=RHOFXSUM(i,j)/WTXSUM(i,j);
            end
        end
    end
    % Vy nodes
    for i=1:1:Ny
        for j=1:1:Nx1
            if(WTYSUM(i,j)>0)
                RHOY(i,j)=RHOYSUM(i,j)/WTYSUM(i,j);
                RHOFY(i,j)=RHOFYSUM(i,j)/WTYSUM(i,j);
            end
        end
    end
    
    % Save viscosity
    if(timestep==1)
        ETA1=ETA0;
        ETA=ETA0;
        ETA50=ETA0;
    else
        ETA1=ETA00;
        ETA=ETA00;
        ETA50=ETA00;
    end
    YNY00=YNY0;
    
    OM=OM0;
    
    % Multiple solving of equations
    if(yndtdecrease==0)
        dt=max(min(dt*dtkoefup,dtelastic0),dtmin);
    else
        dt=max(min(dt,dtelastic0),dtmin);
    end
    yndtdecrease=0;
    dt00=dt;
    DSYLSQ=0;
    ynlast=0;
    for iterstep=1:1:niterglobal
        
        fprintf('iterstep: %d \n',iterstep);
        
        % Limiting viscosity
        etamincur=shearmod*dt*1e-4;
        
        % External P-nodes: symmetry
        pt(:,1)=pt(:,2);
        pt(:,Nx1)=pt(:,Nx);
        pt(1,:)=pt(2,:);
        pt(Ny1,:)=pt(Ny,:);
        pf(:,1)=pf(:,2);
        pf(:,Nx1)=pf(:,Nx);
        pf(1,:)=pf(2,:);
        pf(Ny1,:)=pf(Ny,:);
        % Basic nodes
        for i=1:1:Ny
            for j=1:1:Nx
                if(ETA(i,j)<etamincur)
                    ETA(i,j)=etamincur;
                end
                % Compute plastic strain rate
                if(ETA(i,j)<ETA0(i,j))
                    
                    SXY2=(SXY(i,j)^2);
                    SZX2=(SZX(i,j)^2);
                    SZY2=(SZY(i,j)^2);
                    SIIB=(0.5*(SXX(i,j)^2+SYY(i,j)^2)+SXY2+SZX2+SZY2)^0.5;
                    
                    %SIIB=(SXY(i,j)^2+SZX(i,j)^2+SZY(i,j)^2 ...
                    %    +1/2*((SXX(i,j)+SXX(i+1,j)+SXX(i,j+1)+SXX(i+1,j+1))/4)^2+...
                    %    +1/2*((SYY(i,j)+SYY(i+1,j)+SYY(i,j+1)+SYY(i+1,j+1))/4)^2+...
                    %    +1/2*((-SXX(i,j)-SYY(i,j)-SXX(i+1,j)-SYY(i+1,j)...
                    %    -SXX(i,j+1)-SYY(i,j+1)-SXX(i+1,j+1)-SYY(i+1,j+1))/4)^2)^0.5;
                    
                    EIIB(i,j)=dy/faultwidth*(SIIB/2/ETA(i,j)-SIIB/2/ETA0(i,j));
                    IETAPLB(i,j)=(1/ETA(i,j)-1/ETA0(i,j));
                else
                    EIIB(i,j)=0;
                    IETAPLB(i,j)=0;
                end
            end
        end
        % Computing viscosity and dilatation in pressure nodes
        for i=2:1:Ny
            for j=2:1:Nx
                % Compute viscoplastic viscosity
                IETAPL=(IETAPLB(i-1,j-1)+IETAPLB(i,j-1)+IETAPLB(i-1,j)+IETAPLB(i,j))/4;
                if(YNY0(i-1,j-1)>0 || YNY0(i,j-1)>0 || YNY0(i-1,j)>0 || YNY0(i,j)>0)
                    ETAP(i,j)=1/(1/ETAP0(i,j)+IETAPL);
                    ETAB(i,j)=1/(1/ETAB0(i,j)+dy/faultwidth*IETAPL*POR(i,j));
                else
                    ETAP(i,j)=ETAP0(i,j);
                    ETAB(i,j)=ETAB0(i,j);
                end
                % Check viscosity
                if(ETAP(i,j)<etamincur)
                    ETAP(i,j)=etamincur;
                end
                if(ETAB(i,j)*POR(i,j)<etamincur)
                    ETAB(i,j)=etamincur/POR(i,j);
                end
                % Pores compressibility
                GGGB(i,j)=GGGP(i,j)/POR(i,j);
                % Dilation
                % Zhao and Cai,International Journal of Rock Mechanics & Mining Sciences 47 (2010) 368–384
                % Weak sandstone parameters
                %aa1=20.93;
                %aa2=35.28;
                %aa3=2.34;
                %bb1=0.99;
                %bb2=44.39;
                %bb3=0.73;
                %cc1=0.37;
                %cc2=3.54;
                %cc3=0.47;
                %ss3=min(max((pt(i,j)-pf(i,j))*1e-6,0),100); % SIGMA3, MPa
                %aa=aa1+aa2*exp(-ss3/aa3);
                %bb=bb1+bb2*exp(-ss3/bb3);
                %cc=cc1+cc2/100*ss3^cc3;
                
                %gammapij=plstrain*2*100; % Plastic Strain, %
                %gammapi1j=plstrain*2*100; % Plastic Strain, %
                %gammapij1=plstrain*2*100; % Plastic Strain, %
                %gammapi1j1=plstrain*2*100;% Plastic Strain, %
                %dilij=sin(aa*bb*(exp(-bb*gammapij)-exp(-cc*gammapij))/(cc-bb)/180*pi);
                %dili1j=sin(aa*bb*(exp(-bb*gammapi1j)-exp(-cc*gammapi1j))/(cc-bb)/180*pi);
                %dilij1=sin(aa*bb*(exp(-bb*gammapij1)-exp(-cc*gammapij1))/(cc-bb)/180*pi);
                %dili1j1=sin(aa*bb*(exp(-bb*gammapi1j1)-exp(-cc*gammapi1j1))/(cc-bb)/180*pi);
                DILP(i,j)=0;%2*(dili1j1*EIIB(i-1,j-1)+dilij1*EIIB(i,j-1)+dili1j*EIIB(i-1,j)+dilij*EIIB(i,j))/4;
            end
        end

        gggbkoef=1;
        % Compute pscales
        % ptscale=min(min(ETA))/dx;
        % ptscale=shearmod*dt*1e-3/dx;
        if(dt>2e+5)
            pfscale=min(min(ETAD(2:Ny,1:Nx)))*dx*1e+17/dt^2;
        else
            pfscale=min(min(ETAD(2:Ny,1:Nx)))*dx;
        end
        
        ptscale=pfscale;
        
        % 5)Composing global matrixes L(), R()
        for j=1:1:Nx1
            for i=1:1:Ny1
                % Computing global indexes for vx,vy,p
                kp=((j-1)*Ny1+(i-1))*7+1;
                kx=kp+1;
                ky=kp+2;
                kpf=kp+3;
                kxf=kp+4;
                kyf=kp+5;
                kz=kp+6;
                
                
                % FABIAN (note that kp is with 7 indexes)
                
                % 5o) Composing equation for vzs (out-of-plane component)
                if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
                    
                    % Ghost nodes: 1*vzs=0
                    if(j==Nx1)
                        L(kz,kz)=1;
                        R(kz)=0;
                    end
                    
                    % Upper boundary
                    if(i==1 && j<Nx1)
                        L(kz,kz)=1;
                        L(kz,kz+7)=1;
                        R(kz)=2*bclower;
                    end
                    
                    % Lower boundary
                    if(i==Ny1 && j<Nx1)
                        L(kz,kz)=1;
                        L(kz,kz-7)=1;
                        R(kz)=2*-bclower;
                    end
                    
                    % Left boundary:
                    if(j==1 && i>1 && i<Ny1)
                        L(kz,kz)=1;
                        L(kz,kz+7*Ny1)=-1;
                        R(kz)=0;
                    end
                    
                    % Right boundary
                    if(j==Nx && i>1 && i<Ny1)
                        L(kz,kz)=1;
                        L(kz,kz-7*Ny1)=-1;
                        R(kz)=0;
                    end
                    
                else
                    
                    % Total Z-Stokes: dSIGMAzxt'/dx+dSIGMAzyt'/dy-dPt/dx=-RHOt*gz
                    % SIGMAijt=2*ETA*EPSILONijs*K+SIGMAijt0*(1-K)
                    %             vzs2
                    %              |
                    %        (P)  SZY   (P)
                    %              |
                    %  vzs1--SZX--vzs3--SZX--vzs5
                    %              |
                    %        (P)  SZY   (P)
                    %              |
                    %             vzs4
                    
                    ETAXY=ETA(i,j);
                    % Shear modulus
                    GXY=GGG(i,j);
                    % Viscoelasticity factor
                    KXY=dt*GXY/(dt*GXY+ETAXY);
                    % Numerical viscosity
                    ETAXY=ETAXY*KXY;
                    % Numerical stresses
                    SZY1=SZY0(i-1,j)*(1-KXY);
                    SZY2=SZY0(i,j)*(1-KXY);
                    SZX1=SZX0(i,j-1)*(1-KXY);
                    SZX2=SZX0(i,j)*(1-KXY);
                    % Left part
                    L(kz,kz)=-2*ETAXY/dx^2-2*ETAXY/dy^2-ascale*RHOX(i,j)/dt; %vzs3
                    L(kz,kz-Ny1*7)=ETAXY/dx^2; %vzs1
                    L(kz,kz+Ny1*7)=ETAXY/dx^2; %vzs5
                    L(kz,kz-7)=ETAXY/dy^2; %vzs2
                    L(kz,kz+7)=ETAXY/dy^2; %vzs4
                    % Right part
                    R(kz)=-RHOX(i,j)*(ascale*VZ0(i,j)/dt)-(SZX2-SZX1)/dx-(SZY2-SZY1)/dy;
                end
                
                
                % FABIAN : end!
                
                
                % 5a) Composing equation for vxs
                if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
                    % Ghost nodes: 1*vxs=0
                    if(j==Nx1)
                        L(kx,kx)=1;
                        R(kx)=0;
                    end
                    
                    
                    % Upper boundary
                    if(i==1 && j<Nx1)
                        L(kx,kx)=1;
                        L(kx,kx+7)=0;
                        R(kx)=0;
                    end
                    
                    % Lower boundary
                    if(i==Ny1 && j<Nx1)
                        L(kx,kx)=1;
                        L(kx,kx-7)=0;
                        R(kx)=0;
                    end
                    
                    % Left boundary:
                    if(j==1 && i>1 && i<Ny1)
                        L(kx,kx)=1;
                        L(kx,kx+7*Ny1)=-1;
                        R(kx)=0;
                    end
                    
                    % Right boundary
                    if(j==Nx && i>1 && i<Ny1)
                        L(kx,kx)=1;
                        L(kx,kx-7*Ny1)=-1;
                        R(kx)=0;
                    end
                    
                else
                    
                    % Total X-Stokes: dSIGMAxxt'/dx+dSIGMAxyt'/dy-dPt/dx=-RHOt*gx
                    % SIGMAijt=2*ETA*EPSILONijs*K+SIGMAijt0*(1-K)
                    %             vxs2
                    %        vys1  |    vys3
                    %              |
                    %  vxs1--Pt1--vxs3--Pt2--vxs5
                    %              |
                    %        vys2  |    vys4
                    %             vxs4
                    % Viscosity
                    ETAXY1=ETA(i-1,j);
                    ETAXY2=ETA(i,j);
                    ETAXX1=ETAP(i,j);
                    ETAXX2=ETAP(i,j+1);
                    % Shear modulus
                    GXY1=GGG(i-1,j);
                    GXY2=GGG(i,j);
                    GXX1=GGGP(i,j);
                    GXX2=GGGP(i,j+1);
                    % Viscoelasticity factor
                    KXY1=dt*GXY1/(dt*GXY1+ETAXY1);
                    KXY2=dt*GXY2/(dt*GXY2+ETAXY2);
                    KXX1=dt*GXX1/(dt*GXX1+ETAXX1);
                    KXX2=dt*GXX2/(dt*GXX2+ETAXX2);
                    % Numerical viscosity
                    ETAXY1=ETAXY1*KXY1;
                    ETAXY2=ETAXY2*KXY2;
                    ETAXX1=ETAXX1*KXX1;
                    ETAXX2=ETAXX2*KXX2;
                    % Numerical stresses
                    SXY1=SXY0(i-1,j)*(1-KXY1);
                    SXY2=SXY0(i,j)*(1-KXY2);
                    SXX1=SXX0(i,j)*(1-KXX1);
                    SXX2=SXX0(i,j+1)*(1-KXX2);
                    % Density derivatives
                    dRHOdx=(RHOX(i,j+1)-RHOX(i,j-1))/2/dx;
                    dRHOdy=(RHO(i,j)-RHO(i-1,j))/dy;
                    % Left part
                    L(kx,kx)=-(ETAXX1+ETAXX2)/dx^2 ...
                        -(ETAXY1+ETAXY2)/dy^2-gx*dt*dRHOdx-ascale*RHOX(i,j)/dt; %vxs3                    
                    L(kx,kx-Ny1*7)=ETAXX1/dx^2; %vxs1
                    L(kx,kx+Ny1*7)=ETAXX2/dx^2; %vxs5
                    L(kx,kx-7)=ETAXY1/dy^2; %vxs2
                    L(kx,kx+7)=ETAXY2/dy^2; %vxs4
                    L(kx,ky-7)=ETAXY1/dx/dy-ETAXX1/dx/dy-gx*dt*dRHOdy/4; %vys1
                    L(kx,ky)=-ETAXY2/dx/dy+ETAXX1/dx/dy-gx*dt*dRHOdy/4; %vys2
                    L(kx,ky-7+Ny1*7)=-ETAXY1/dx/dy+ETAXX2/dx/dy-gx*dt*dRHOdy/4; %vys3
                    L(kx,ky+Ny1*7)=ETAXY2/dx/dy-ETAXX2/dx/dy-gx*dt*dRHOdy/4; %vys4
                    L(kx,kp)=ptscale/dx; %Pt1'
                    L(kx,kp+Ny1*7)=-ptscale/dx; %Pt2'
                    % Right part
                    R(kx)=-RHOX(i,j)*(ascale*VX0(i,j)/dt+gx)-(SXX2-SXX1)/dx-(SXY2-SXY1)/dy;
                end
                
                % 5b) Composing equation for vys
                if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
                    % Ghost nodes: 1*vys=0
                    if(i==Ny1)
                        L(ky,ky)=1;
                        R(ky)=0;
                    end
                    
                    % Left boundary
                    % Free Slip
                    if(j==1)
                        L(ky,ky)=1;
                        L(ky,ky+Ny1*7)=1;
                        R(ky)=0;
                    end
                    
                    % Right boundary
                    % Free Slip
                    if(j==Nx1)
                        L(ky,ky)=1;
                        L(ky,ky-Ny1*7)=1;
                        R(ky)=0;
                    end
                    
                    % Upper boundary: no penetration
                    if(i==1 && j>1 && j<Nx1)
                        L(ky,ky)=1;
                        R(ky)=0;
                    end
                    
                    % Lower boundary: no penetration
                    if(i==Ny && j>1 && j<Nx1)
                        L(ky,ky)=1;
                        R(ky)=0;
                    end
                    
                else
                    % Total Y-Stokes: dSIGMAyxt'/dx+dSIGMAyyt'/dy-dPt/dy=-RHOt*gy
                    % y-Stokes equation: dSIGMA'yx/dx+dSIGMA'yy/dy-dP/dy=-RHO*gy
                    %
                    %               vys2
                    %                |
                    %         vxs1  Pt1  vxs3
                    %                |
                    %   vys1---------vys3--------vys5
                    %                |
                    %         vxs2  Pt2  vxs4
                    %                |
                    %               vys4
                    % Viscosity
                    ETAXY1=ETA(i,j-1);
                    ETAXY2=ETA(i,j);
                    ETAYY1=ETAP(i,j);
                    ETAYY2=ETAP(i+1,j);
                    % Shear modulus
                    GXY1=GGG(i,j-1);
                    GXY2=GGG(i,j);
                    GYY1=GGGP(i,j);
                    GYY2=GGGP(i+1,j);
                    % Viscoelasticity factor
                    KXY1=dt*GXY1/(dt*GXY1+ETAXY1);
                    KXY2=dt*GXY2/(dt*GXY2+ETAXY2);
                    KYY1=dt*GYY1/(dt*GYY1+ETAYY1);
                    KYY2=dt*GYY2/(dt*GYY2+ETAYY2);
                    % Numerical viscosity
                    ETAXY1=ETAXY1*KXY1;
                    ETAXY2=ETAXY2*KXY2;
                    ETAYY1=ETAYY1*KYY1;
                    ETAYY2=ETAYY2*KYY2;
                    % Numerical stresses
                    SXY1=SXY0(i,j-1)*(1-KXY1);
                    SXY2=SXY0(i,j)*(1-KXY2);
                    SYY1=SYY0(i,j)*(1-KYY1);
                    SYY2=SYY0(i+1,j)*(1-KYY2);
                    % Density derivatives
                    dRHOdy=(RHOY(i+1,j)-RHOY(i-1,j))/2/dy;
                    dRHOdx=(RHO(i,j)-RHO(i,j-1))/dx;
                    % Left part
                    L(ky,ky)=-(ETAYY1+ETAYY2)/dy^2-...
                        (ETAXY1+ETAXY2)/dx^2-gy*dt*dRHOdy-ascale*RHOY(i,j)/dt; %vys3
                    L(ky,ky-Ny1*7)=ETAXY1/dx^2; %vys1
                    L(ky,ky+Ny1*7)=ETAXY2/dx^2; %vys5
                    L(ky,ky-7)=ETAYY1/dy^2; %vys2
                    L(ky,ky+7)=ETAYY2/dy^2; %vys4
                    L(ky,kx-Ny1*7)=ETAXY1/dx/dy-ETAYY1/dx/dy-gy*dt*dRHOdx/4; %vxs1
                    L(ky,kx+7-Ny1*7)=-ETAXY1/dx/dy+ETAYY2/dx/dy-gy*dt*dRHOdx/4; %vxs2
                    L(ky,kx)=-ETAXY2/dx/dy+ETAYY1/dx/dy-gy*dt*dRHOdx/4; %vxs3
                    L(ky,kx+7)=ETAXY2/dx/dy-ETAYY2/dx/dy-gy*dt*dRHOdx/4; %vxs4
                    L(ky,kp)=ptscale/dy; %Pt1'
                    L(ky,kp+7)=-ptscale/dy; %Pt2'
                    % Right part
                    R(ky)=-RHOY(i,j)*(ascale*VY0(i,j)/dt+gy)-(SYY2-SYY1)/dy-(SXY2-SXY1)/dx;
                end
                
                
                % 5c) Composing equation for Pt
                if(i==1 || j==1 || i==Ny1 || j==Nx1)
                    % BC equation: 1*Pt=0
                    L(kp,kp)=1;
                    R(kp)=0;
                    
                else
                    % Solid Continuity: dVxs/dx+dVys/dy+(Pt-Pf)/ETAbulk=0
                    %              vys1
                    %               |
                    %        vxs1--Pt,Pf--vxs2
                    %               |
                    %              vys2
                    % Drained compressibility
                    BETTADRAINED=(1/GGGB(i,j)+BETTASOLID)/(1-POR(i,j));
                    % Biott-Willis koefficient
                    KBW=1-BETTASOLID/BETTADRAINED;
                    % Left part
                    L(kp,kx-Ny1*7)=-1/dx; %vxs1
                    L(kp,kx)=1/dx; %vxs2
                    L(kp,ky-7)=-1/dy; %vys1
                    L(kp,ky)=1/dy; %vys2
                    L(kp,kp)=ptscale*(1/ETAB(i,j)/(1-POR(i,j))+gggbkoef*BETTADRAINED/dt); %Pt
                    L(kp,kpf)=-pfscale*(1/ETAB(i,j)/(1-POR(i,j))+gggbkoef*BETTADRAINED*KBW/dt); %Pf
                    % Right part
                    R(kp)=gggbkoef*BETTADRAINED*(PT0(i,j)-KBW*PF0(i,j))/dt+DILP(i,j);
                end
                
                % 5d) Composing equation for vxD
                if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
                    % Ghost nodes: 1*vxs=0
                    if(j==Nx1)
                        L(kxf,kxf)=1;
                        R(kxf)=0;
                    end
                    
                    % Upper boundary: symmetry
                    if(i==1 && j<Nx1)
                        L(kxf,kxf)=1;
                        L(kxf,kxf+7)=-1;
                        R(kxf)=0;
                    end
                    
                    % Lower boundary: symmetry
                    if(i==Ny1 && j<Nx1)
                        L(kxf,kxf)=1;
                        L(kxf,kxf-7)=-1;
                        R(kxf)=0;
                    end
                    
                    % Left boundary
                    % no penetration
                    if(j==1)
                        L(kxf,kxf)=1;
                        R(kxf)=0;%bcvxfleft;
                    end
                    
                    % Right boundary
                    % no penetration
                    if(j==Nx)
                        L(kxf,kxf)=1;
                        R(kxf)=0;
                    end
                    
                else
                    % Fluid X-Darsi: -ETAfluid/K*VxD-dPf/dx=-RHOf*gx+RHOf*DVxs/Dt
                    %
                    %  Pf1---vxD,vxs---Pf2
                    %
                    % Left part
                    L(kxf,kxf)=-ETAD(i,j)-RHOFX(i,j)/POR(i,j)*ascale/dt; %vxD
                    L(kxf,kx)=-RHOFX(i,j)*ascale/dt; %vxs
                    L(kxf,kpf)=pfscale/dx; %Pf1'
                    L(kxf,kpf+Ny1*7)=-pfscale/dx; %Pf2'
                    % Right part
                    R(kxf)=-RHOFX(i,j)*(ascale*VXF0(i,j)/dt+gx);
                end
                
                % 5e) Composing equation for vyD
                if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
                    % Ghost nodes: 1*vxs=0
                    if(i==Ny1)
                        L(kyf,kyf)=1;
                        R(kyf)=0;
                    end
                    
                    % Left boundary
                    % symmetry
                    if(j==1 && i>1 && i<Ny)
                        L(kyf,kyf)=1;
                        L(kyf,kyf+Ny1*7)=-1;
                        R(kyf)=0;
                    end
                    
                    % Right boundary
                    % symmetry
                    if(j==Nx1 && i>1 && i<Ny)
                        L(kyf,kyf)=1;
                        L(kyf,kyf-Ny1*7)=-1;
                        R(kyf)=0;
                    end
                    
                    % Upper boundary: no penetration
                    if(i==1)
                        L(kyf,kyf)=1;
                        R(kyf)=bcvyflower;
                    end
                    
                    % Lower boundary: no penetration
                    if(i==Ny)
                        L(kyf,kyf)=1;
                        R(kyf)=bcvyflower;
                    end
                else
                    % Fluid Y-Darsi: -ETAfluid/K*VyD-dPf/dy=-RHOf*gy+RHOf*DVys/Dt
                    %
                    %   Pf1
                    %    |
                    %   vyD,vy
                    %    |
                    %   Pf2                |
                    %
                    % Left part
                    L(kyf,kyf)=-ETAD(i,j)-RHOFY(i,j)/POR(i,j)*ascale/dt; %vyD
                    L(kyf,ky)=-RHOFY(i,j)*ascale/dt; %vys
                    L(kyf,kpf)=pfscale/dy; %Pf1'
                    L(kyf,kpf+7)=-pfscale/dy; %Pf2'
                    % Right part
                    R(kyf)=-RHOFY(i,j)*(ascale*VYF0(i,j)/dt+gy);
                end
                
                
                % 5f) Composing equation for Pf
                if(i==1 || j==1 || i==Ny1 || j==Nx1 || i==2 || i==Ny)
                    % BC equation: 1*Pf=0
                    L(kpf,kpf)=1;
                    R(kpf)=0;
                    
                    % Real BC
                    if(i==2 && (j>2 || gggbkoef==1))
                        L(kpf,kpf)=1*pfscale;
                        L(kpf,kp)=-1*ptscale;
                        R(kpf)=-PTFDIFF(i,j);
                    end
                    
                    % Real BC
                    if(i==Ny && (j<Nx || gggbkoef==1))
                        L(kpf,kpf)=1*pfscale;
                        L(kpf,kp)=-1*ptscale;
                        R(kpf)=-PTFDIFF(i,j);
                    end

                    
                else
                    % Fluid Continuity: dVxD/dx+dVyD/dy-(Pt-Pf)/ETAbulk=0
                    %              vyD1
                    %               |
                    %        vxD1--Pt,Pf--vxD2
                    %               |
                    %              vyD2
                    % Compute elastic coefficients
                    % Drained compressibility
                    BETTADRAINED=(1/GGGB(i,j)+BETTASOLID)/(1-POR(i,j));
                    % Biott-Willis koefficient
                    KBW=1-BETTASOLID/BETTADRAINED;
                    % Skempton koefficient
                    KSK=(BETTADRAINED-BETTASOLID)/(BETTADRAINED-BETTASOLID+POR(i,j)*(BETTAFLUID-BETTASOLID));
                    % Left part
                    L(kpf,kxf-Ny1*7)=-1/dx; %vxs1
                    L(kpf,kxf)=1/dx; %vxs2
                    L(kpf,kyf-7)=-1/dy; %vys1
                    L(kpf,kyf)=1/dy; %vys2
                    L(kpf,kp)=-ptscale*(1/ETAB(i,j)/(1-POR(i,j))+gggbkoef*BETTADRAINED*KBW/dt); %Pt
                    L(kpf,kpf)=pfscale*(1/ETAB(i,j)/(1-POR(i,j))+gggbkoef*BETTADRAINED*KBW/KSK/dt); %Pf
                    % Right part
                    R(kpf)=-gggbkoef*BETTADRAINED*KBW*(PT0(i,j)-1/KSK*PF0(i,j))/dt-DILP(i,j);
                end
            end
        end
        
        % 6) Solving matrix
        S=L\R;
        
        
        % 7) Reload solution
        for j=1:1:Nx1
            for i=1:1:Ny1
                % Global indexes for vx,vy,P
                kp=((j-1)*Ny1+(i-1))*7+1;
                kx=kp+1;
                ky=kp+2;
                kpf=kp+3;
                kxf=kp+4;
                kyf=kp+5;
                kz=kp+6;
                % Reload solution
                pt(i,j)=S(kp)*ptscale;
                vxs(i,j)=S(kx);
                vys(i,j)=S(ky);
                vzs(i,j)=S(kz);
                pf(i,j)=S(kpf)*pfscale;
                vxD(i,j)=S(kxf);
                vyD(i,j)=S(kyf);
            end
        end

        Vmax = max(max(VSLIPB));

        % Velocity change
        % FABIAN
        vzs(:,Nx1)=vzs(:,Nx);
        DVX0=vxs-VX0;
        DVY0=vys-VY0;
        
        % FABIAN
        DVZ0=vzs-VZ0;
        
        % Define timestep
        dt0=dt;
        yn=0;
        
        % Plastic iterations
        % Compute strain rate, stress and stress change
        ESP=zeros(Ny,Nx);
        
        EXY=zeros(Ny,Nx);
        EZX=zeros(Ny,Nx);
        EZY=zeros(Ny,Nx);
        
        SXY=zeros(Ny,Nx);
        
        % FABIAN
        SZX=zeros(Ny,Nx);
        SZY=zeros(Ny,Nx);
                
        %DSXY=zeros(Ny,Nx);
        EXX=zeros(Ny1,Nx1);
        SXX=zeros(Ny1,Nx1);
        DSXX=zeros(Ny1,Nx1);
        EYY=zeros(Ny1,Nx1);
        SYY=zeros(Ny1,Nx1);
        DSYY=zeros(Ny1,Nx1);
        EII=zeros(Ny1,Nx1);
        %EIIVP=zeros(Ny1,Nx1);
        SII=zeros(Ny1,Nx1);
        %DIS=zeros(Ny1,Nx1);
        
        EL_DECOM=zeros(Ny1,Nx1);   % Elastic (de)compaction
        VIS_COMP=zeros(Ny1,Nx1);   % Viscous compaction
        
        % Process internal basic nodes
        for i=1:1:Ny
            for j=1:1:Nx
                
                ESP(i,j)=1/2*((vys(i,j+1)-vys(i,j))/dx-(vxs(i+1,j)-vxs(i,j))/dy);
                EXY(i,j)=1/2*((vxs(i+1,j)-vxs(i,j))/dy+(vys(i,j+1)-vys(i,j))/dx);
                
                % FABIAN
                EZX(i,j)=1/2*((vzs(i,j+1)-vzs(i,j))/dx+(vxs(i+1,j)-vxs(i,j))/dx);
                EZY(i,j)=1/2*((vzs(i+1,j)-vzs(i,j))/dy+(vys(i,j+1)-vys(i,j))/dy);
                
                KXY=dt*GGG(i,j)/(dt*GGG(i,j)+ETA(i,j));

                
                SXY(i,j)=2*ETA(i,j)*EXY(i,j)*KXY+SXY0(i,j)*(1-KXY);
                
                % FABIAN
                SZX(i,j)=2*ETA(i,j)*EZX(i,j)*KXY+SZX0(i,j)*(1-KXY);
                SZY(i,j)=2*ETA(i,j)*EZY(i,j)*KXY+SZY0(i,j)*(1-KXY);
            end
        end
                
        % Process pressure cells
        for i=2:1:Ny
            for j=2:1:Nx
                % EXX, SXX, DSXX
                EXX(i,j)=(2*(vxs(i,j)-vxs(i,j-1))/dx-(vys(i,j)-vys(i-1,j))/dy)/3;
                EYY(i,j)=(2*(vys(i,j)-vys(i-1,j))/dy-(vxs(i,j)-vxs(i,j-1))/dx)/3;
                KXX=dt*GGGP(i,j)/(dt*GGGP(i,j)+ETAP(i,j));
                                
                %SXX(i,j)=2*GGG(i,j)*dt*EXX(i,j)+SXX0(i,j);
                %SYY(i,j)=2*GGG(i,j)*dt*EYY(i,j)+SYY0(i,j);
                
                SXX(i,j)=2*ETAP(i,j)*EXX(i,j)*KXX+SXX0(i,j)*(1-KXX);
                SYY(i,j)=2*ETAP(i,j)*EYY(i,j)*KXX+SYY0(i,j)*(1-KXX);
                
                DSXX(i,j)=SXX(i,j)-SXX0(i,j);
                DSYY(i,j)=SYY(i,j)-SYY0(i,j);
            end
        end
       
        % External P-nodes: symmetry
        pt(:,1)=pt(:,2);
        pt(:,Nx1)=pt(:,Nx);
        pt(1,:)=pt(2,:);
        pt(Ny1,:)=pt(Ny,:);
        pf(:,1)=pf(:,2);
        pf(:,Nx1)=pf(:,Nx);
        pf(1,:)=pf(2,:);
        pf(Ny1,:)=pf(Ny,:);
        EXX(:,1)=EXX(:,2);
        EXX(:,Nx1)=EXX(:,Nx);
        EXX(1,:)=EXX(2,:);
        EXX(Ny1,:)=EXX(Ny,:);
        SXX(:,1)=SXX(:,2);
        SXX(:,Nx1)=SXX(:,Nx);
        SXX(1,:)=SXX(2,:);
        SXX(Ny1,:)=SXX(Ny,:);
        SXX0(:,1)=SXX0(:,2);
        SXX0(:,Nx1)=SXX0(:,Nx);
        SXX0(1,:)=SXX0(2,:);
        SXX0(Ny1,:)=SXX0(Ny,:);
        EYY(:,1)=EYY(:,2);
        EYY(:,Nx1)=EYY(:,Nx);
        EYY(1,:)=EYY(2,:);
        EYY(Ny1,:)=EYY(Ny,:);
        SYY(:,1)=SYY(:,2);
        SYY(:,Nx1)=SYY(:,Nx);
        SYY(1,:)=SYY(2,:);
        SYY(Ny1,:)=SYY(Ny,:);
        
        % FABIAN
        SZX(:,Nx)=SZX(:,Nx-1);
        EZX(:,Nx)=EZX(:,Nx-1);
        
        SYY0(:,1)=SYY0(:,2);
        SYY0(:,Nx1)=SYY0(:,Nx);
        SYY0(1,:)=SYY0(2,:);
        SYY0(Ny1,:)=SYY0(Ny,:);
        ETAP(:,1)=ETAP(:,2);
        ETAP(:,Nx1)=ETAP(:,Nx);
        ETAP(1,:)=ETAP(2,:);
        ETAP(Ny1,:)=ETAP(Ny,:);
        ETAB(:,1)=ETAB(:,2);
        ETAB(:,Nx1)=ETAB(:,Nx);
        ETAB(1,:)=ETAB(2,:);
        ETAB(Ny1,:)=ETAB(Ny,:);
        GGGP(:,1)=GGGP(:,2);
        GGGP(:,Nx1)=GGGP(:,Nx);
        GGGP(1,:)=GGGP(2,:);
        GGGP(Ny1,:)=GGGP(Ny,:);
        GGGB(:,1)=GGGB(:,2);
        GGGB(:,Nx1)=GGGB(:,Nx);
        GGGB(1,:)=GGGB(2,:);
        GGGB(Ny1,:)=GGGB(Ny,:);
                
        % Compute stress and strain rate invariants and dissipation
        % Process pressure cells
        for i=2:1:Ny
            for j=2:1:Nx
                
                % EXY term is averaged from four surrounding basic nodes
                EXY2=(EXY(i,j)^2+EXY(i-1,j)^2+EXY(i,j-1)^2+EXY(i-1,j-1)^2)/4;
                
                % FABIAN
                EZX2=(EZX(i,j)^2+EZX(i-1,j)^2+EZX(i,j-1)^2+EZX(i-1,j-1)^2)/4;
                EZY2=(EZY(i,j)^2+EZY(i-1,j)^2+EZY(i,j-1)^2+EZY(i-1,j-1)^2)/4;
                EII(i,j)=(0.5*(EXX(i,j)^2+EYY(i,j)^2)+EXY2+EZX2+EZY2)^0.5;
                
                
                % Second strain rate invariant SII
                % SXY term is averaged from four surrounding basic nodes
                SXY2=(SXY(i,j)^2+SXY(i-1,j)^2+SXY(i,j-1)^2+SXY(i-1,j-1)^2)/4;
                
                % FABIAN
                SZX2=(SZX(i,j)^2+SZX(i-1,j)^2+SZX(i,j-1)^2+SZX(i-1,j-1)^2)/4;
                SZY2=(SZY(i,j)^2+SZY(i-1,j)^2+SZY(i,j-1)^2+SZY(i-1,j-1)^2)/4;
                SII(i,j)=(0.5*(SXX(i,j)^2+SYY(i,j)^2)+SXY2+SZX2+SZY2)^0.5;
                
            end
        end
        
        % Update viscosity for yielding
        AXY=zeros(Ny,Nx);
        % dt0=dt;dt=dt*1.1;
        ETA5=ETA0;
        % Basic nodes
        DSY=zeros(Ny,Nx);
        YNY=zeros(Ny,Nx);
        SigmaY=zeros(Ny,Nx);
        VSLIPB=zeros(Ny,Nx);
        SII_fault=zeros(Ny,Nx);
        ynpl=0;
        ddd=0;
        dtrobert=dt;
        dtlapusta=3e7;
        OM5=OM;
        
        % Power law plasticity model Yi et al., 2018
        % Journal of Offshore Mechanics and Arctic Engineering
        dtslip=1e+30;
        if(timestep>tyield)
            for i=1:1:Ny
                for j=1:1:Nx                    
                    %if(y(i)>=0) %y(line_fault)-(dy) && y(i)<=y(line_fault)+(dy))
                    if(i>=line_fault-1 && i<=line_fault+1)          
	            %if i==line_fault

                        SXY2=(SXY(i,j)^2);
                        % FABIAN
                        SZX2=(SZX(i,j)^2);
                        SZY2=(SZY(i,j)^2);
                        SIIB(i,j)=(0.5*(SXX(i,j)^2+SYY(i,j)^2)+SXY2+SZX2+SZY2)^0.5;                        
                     
                        ptB=pt(i,j);%(pt(i,j)+pt(i+1,j)+pt(i,j+1)+pt(i+1,j+1))/4;
                        pfB=pf(i,j);%(pf(i,j)+pf(i+1,j)+pf(i,j+1)+pf(i+1,j+1))/4;	
                        % Computing "elastic" stress invariant
                        kfxy=ETA(i,j)/(GGG(i,j)*dt+ETA(i,j));
                        siiel=SIIB(i,j)/kfxy;
                        % Compute maximal stress invariant with viscous viscosity
                        kfxy0=ETA0(i,j)/(GGG(i,j)*dt+ETA0(i,j));
                        SIIB0=siiel*kfxy0;
                        
                        % Assign PEFFB
                        PEFFB0(i,j)=(ptB-pfB);
                        PEFFB(i,j)=(ptB-pfB);
                        PEFFB1(i,j)=PEFFB(i,j)-PEFFB0(i,j);
                        % Compute old viscoplastic slip rate
                        % Compute PEFF
                        prB=(ptB-pfB);
                        
                        if(prB<1e6)
                            prB=1e6;
                        end
                        % Compute old power law strain rate
                        SIIB1=SIIB(i,j);
                        
                        %Compute slip velocity for current stress invariant and state
                        V=2*V0*sinh(max((SIIB1),0)/ARSF(i,j)/prB)*...
                            exp(-(BRSF(i,j)*OM(i,j)+FRIC(i,j))/ARSF(i,j));
                        
                        EIISLIP=V/dx/2;
                        
                        % Compute new ETAVP
                        ETAPL=SIIB1/2/EIISLIP;
                        ETAVP=1/(1/ETA0(i,j)+1/ETAPL);
                        % Compute new stress invariant
                        kfxy1=ETAVP/(GGG(i,j)*dt+ETAVP);
                        SIIB2=siiel*kfxy1;
                        DSIIB1=SIIB2-SIIB1;
                        
                        
                        %Compute slip velocity for current stress invariant and state
                        V=2*V0*sinh(max((SIIB2),0)/ARSF(i,j)/prB)*...
                            exp(-(BRSF(i,j)*OM(i,j)+FRIC(i,j))/ARSF(i,j));
                        
                        EIISLIP=V/dx/2;
                        
                        % Compute new ETAVP
                        ETAPL=SIIB2/2/EIISLIP;
                        ETAVP=1/(1/ETA0(i,j)+1/ETAPL);
                        % Compute new stress invariant
                        kfxy1=ETAVP/(GGG(i,j)*dt+ETAVP);
                        SIIB3=siiel*kfxy1;
                        DSIIB2=SIIB3-SIIB2;
                        
                        if((DSIIB1>=0 && DSIIB2<=0) || (DSIIB1<=0 && DSIIB2>=0))
                            DSIIB=1e+9;
                            ijk=0;
                            while(abs(DSIIB)>1e-3)
                                SIIB4=(SIIB1+SIIB2)/2;
                                
                                %Compute slip velocity for current stress invariant and state
                                V=2*V0*sinh(max((SIIB4),0)/ARSF(i,j)/prB)*...
                                    exp(-(BRSF(i,j)*OM(i,j)+FRIC(i,j))/ARSF(i,j));
                                
                                EIISLIP=V/dx/2;
                                
                                % Compute new ETAVP
                                ETAPL=SIIB4/2/EIISLIP;
                                ETAVP=1/(1/ETA0(i,j)+1/ETAPL);
                                % Compute new stress invariant
                                kfxy1=ETAVP/(GGG(i,j)*dt+ETAVP);
                                SIIB5=siiel*kfxy1;
                                DSIIB=SIIB5-SIIB4;
                                if((DSIIB>=0 && DSIIB1>=0) || (DSIIB<=0 && DSIIB1<=0))
                                    SIIB1=SIIB4;
                                else
                                    SIIB2=SIIB4;
                                end
                                ijk=ijk+1;
                            end
                        end
                        
                        
                        if(V*dt/LRSF(i,j)>1e-6)
                            OM5(i,j)=log(V0/V+(exp(OM0(i,j))-V0/V)*exp(-V*dt/LRSF(i,j)));
                        else
                            OM5(i,j)=log(exp(OM0(i,j))*(1-V*dt/LRSF(i,j))+V0*dt/LRSF(i,j));
                        end
                        
                        
                        kfxy1=ETAVP/(GGG(i,j)*dt+ETAVP);
                        SIGMA2=siiel*kfxy1;
                        
                        % Compute yielding stress
                        syield=max(syieldmin,(ptB-pfB)*ARSF(i,j)*asinh(V/2/V0*exp((BRSF(i,j)*OM5(i,j)+FRIC(i,j))/ARSF(i,j))));
                        
                        % Compute visco-plastic viscosity
                        D=1;
                        etapl=ETA0(i,j)*syield/(ETA0(i,j)*V/D+syield);
                        
                        %save syield
                        SigmaY(i,j)=syield;
                        VSLIPB(i,j)=V;
                        SII_fault(i,j)=SIIB4;
                        
                        
                        % Timestep criterion, Lapusta et al., 2000; Lapusta and Liu, 2009
                        B=1/BETTASOLID;
                        vi=(3*B-2*GGG(i,j))/(6*B+2*GGG(i,j));
                        k=2/pi*GGG(i,j)/(1-vi)/dx;
                        xi=1/4*(k*LRSF(i,j)/ARSF(i,j)/prB-(BRSF(i,j)-...
                            ARSF(i,j))/ARSF(i,j))^2-k*LRSF(i,j)/ARSF(i,j)/prB;
                        if(xi<0)
                            dTETAmax=min(1-(BRSF(i,j)-ARSF(i,j))*prB/(k*LRSF(i,j)),0.2);
                        else
                            dTETAmax=min(ARSF(i,j)*prB/(k*LRSF(i,j)-(BRSF(i,j)-ARSF(i,j))*prB),0.2);
                        end
                        dtlapusta=min(dtlapusta,dTETAmax*LRSF(i,j)/V);
                        
                        
                        A=syield/siiel;
                        AXY(i,j)=A;
                        % Count old yelding nodes
                        ynn=0;
                        if(YNY0(i,j)>0)
                            ynn=1;
                            DSY(i,j)=SIIB(i,j)-syield;
                            ddd=ddd+DSY(i,j)^2;
                            ynpl=ynpl+1;
                        end
                        % Update viscosity
                        if(A<1)
                            % New viscosity for the basic node
                            etapl=dt*GGG(i,j)*A/(1-A);
                            if(etapl<ETA0(i,j))
                                % Update plastic nodes
                                ETA5(i,j)=etapl^(1-etawt)*ETA(i,j)^etawt;
                                YNY(i,j)=1;
                                % Count yelding nodes
                                if(ynn==0)
                                    DSY(i,j)=SIIB(i,j)-syield;
                                    ddd=ddd+DSY(i,j)^2;
                                    ynpl=ynpl+1;
                                end
                            else
                                ETA5(i,j)=ETA0(i,j);
                            end
                        else
                            ETA5(i,j)=ETA0(i,j);
                        end
                    end
                end
            end
        end
        
        % Compute Error
        DSYLSQ(iterstep)=0;if(ynpl>0);DSYLSQ(iterstep)=(ddd/ynpl)^0.5;end
        if(ynpl==0)
            ETA=ETA0;
        end
        
        % Adjust timestep
        dtpl=dt;
        if (ynpl>0 && iterstep<niterglobal && ynlast>=dtstep && (ynlast>ynlastmax || log10(DSYLSQ(iterstep)/DSYLSQ(iterstep-1))>=0 || log10(DSYLSQ(iterstep)/DSYLSQ(iterstep-1))>log10(errmin/DSYLSQ(iterstep))/(ynlastmax-ynlast)))
            dtpl=dt/dtkoef;
            yn=1;
        end
        
        % Define displacement timesteps
        if(iterstep>1)
            maxvxyz0=maxvxyz;
        end
        maxvxs=max(max(abs(vxs)));
        maxvys=max(max(abs(vys)));
        % FABIAN
        maxvzs=max(max(abs(vzs)));
        maxvxyz=((max(max(vxs))-min(min(vxs)))^2+(max(max(vys))-min(min(vys)))^2+(max(max(vzs))-min(min(vzs)))^2)^0.5;
        stpmaxcur=stpmax1;
        dtx=dt;
        if(dt>dx*stpmaxcur/maxvxs)
            dtx=dx/dtkoefv*stpmaxcur/maxvxs;
            yn=1;
        end
        dty=dt;
        if(dt>dy*stpmaxcur/maxvys)
            dty=dy/dtkoefv*stpmaxcur/maxvys;
            yn=1;
        end
        % FABIAN
        dtz=dt;
        if(dt>dy*stpmaxcur/maxvzs)
            dtz=dx/dtkoefv*stpmaxcur/maxvzs;
            yn=1;
        end
        
        maxvxs=0;
        maxvys=0;
        maxvzs=0;

        for i=1:1:Ny+1
            for j=1:1:Nx+1
                if(i>=bulk_up && i<=bulk_up)
                    maxvxs=max(maxvxs,abs(vxs(i,j)));
                    maxvys=max(maxvys,abs(vys(i,j)));
                    % FABIAN
                    maxvzs=max(maxvzs,abs(vzs(i,j)));
                end
            end
        end
        
        stpmaxcur=stpmax;
        if(dt>dx*stpmaxcur/maxvxs)
            dtx=dx/dtkoefv*stpmaxcur/maxvxs;
            yn=1;
        end
        if(dt>dy*stpmaxcur/maxvys)
            dty=dy/dtkoefv*stpmaxcur/maxvys;
            yn=1;
        end
        % FABIAN
        if(dt>dy*stpmaxcur/maxvzs)
            dtz=dy/dtkoefv*stpmaxcur/maxvzs;
            yn=1;
        end
        
        dtslip=1e30;
        for i=1:1:Ny
            for j=1:1:Nx
                if(VSLIPB(i,j)>0)
                    dtslip=min(dtslip,dx*stpmax/VSLIPB(i,j));
                end
            end
        end
        
        if(ynpl>0 && dtslip<dt)
            yn=1;
            dtslip=dtslip/dtkoefv;
        end
        
        
        % Chose minimal timestep
        if(yn>0 && dt>dtmin)
            dtold=dt;
            dt=max(min([dtx dty dtz dtpl dtslip dtlapusta]),dtmin);
            if(dt<dtold)
                ynlast=0;
            end
        else
            yn=0;
        end
        
        % Exit iterations
        ynstop=0;
        % Velocity change ratio
        if(iterstep>1)
            % FABIAN
            vratio=log10(maxvxyz/maxvxyz0);
        end
        %if(yn==0 && (ynpl>0 && (DSYLSQ(iterstep)<errmin && iterstep>1 && abs(vratio)<vratiomax)))
        if(yn==0 && ((DSYLSQ(iterstep)<errmin && iterstep>1 && abs(vratio)<vratiomax)))
	    ynstop=1;
        else
            % Recomputing ETA
            for i=1:1:Ny
                for j=1:1:Nx
                    ETA(i,j)=max(min(ETA5(i,j),ETA0(i,j)),etamin);
                end
            end
            % Save current viscosity
            ETA50=ETA;
            YNY0=YNY;
            OM=OM5;
        end
        
        % Exit iteration
        if(ynstop==1)
            break;
        end
        ynlast=ynlast+1;
    end
    % Mark dt decrease
    if(dt00>dt)
        yndtdecrease=1;
    end
    
    % Save current viscosity
    ETA00=ETA50;
    OM0=OM;
    % Recheck displacement timestep
    dt=min([dtx dty dtz dt dtlapusta]);
    
    
    % Compute strain rate, stress and stress change
    ESP=zeros(Ny,Nx);
    EXY=zeros(Ny,Nx);
    SXY=zeros(Ny,Nx);
    DSXY=zeros(Ny,Nx);
    EXX=zeros(Ny1,Nx1);
    SXX=zeros(Ny1,Nx1);
    DSXX=zeros(Ny1,Nx1);
    EYY=zeros(Ny1,Nx1);
    SYY=zeros(Ny1,Nx1);
    DSYY=zeros(Ny1,Nx1);
    EII=zeros(Ny1,Nx1);
    EIIVP=zeros(Ny1,Nx1);
    SII=zeros(Ny1,Nx1);
    DSII=zeros(Ny1,Nx1);
    DIS=zeros(Ny1,Nx1);
    DISXY=zeros(Ny1,Nx1);
    DISZX=zeros(Ny1,Nx1);
    DISZY=zeros(Ny1,Nx1);
    
    % Process internal basic nodes
    for i=1:1:Ny
        for j=1:1:Nx
            ESP(i,j)=1/2*((vys(i,j+1)-vys(i,j))/dx-(vxs(i+1,j)-vxs(i,j))/dy);
            EXY(i,j)=1/2*((vxs(i+1,j)-vxs(i,j))/dy+(vys(i,j+1)-vys(i,j))/dx);
            KXY=dt*GGG(i,j)/(dt*GGG(i,j)+ETA(i,j));
            
            % FABIAN
            EZX(i,j)=1/2*((vzs(i,j+1)-vzs(i,j))/dx+(vxs(i+1,j)-vxs(i,j))/dx);
            EZY(i,j)=1/2*((vzs(i+1,j)-vzs(i,j))/dy+(vys(i,j+1)-vys(i,j))/dy);
            
            SXY(i,j)=2*ETA(i,j)*EXY(i,j)*KXY+SXY0(i,j)*(1-KXY);
            
            % FABIAN
            SZX(i,j)=2*ETA(i,j)*EZX(i,j)*KXY+SZX0(i,j)*(1-KXY);
            SZY(i,j)=2*ETA(i,j)*EZY(i,j)*KXY+SZY0(i,j)*(1-KXY);
            
            DSXY(i,j)=SXY(i,j)-SXY0(i,j);
        end
    end
        
    % Process pressure cells
    for i=2:1:Ny
        for j=2:1:Nx
            % EXX, SXX, DSXX
            EXX(i,j)=(2*(vxs(i,j)-vxs(i,j-1))/dx-(vys(i,j)-vys(i-1,j))/dy)/3;
            EYY(i,j)=(2*(vys(i,j)-vys(i-1,j))/dy-(vxs(i,j)-vxs(i,j-1))/dx)/3;
            KXX=dt*GGGP(i,j)/(dt*GGGP(i,j)+ETAP(i,j));
            
            SXX(i,j)=2*ETAP(i,j)*EXX(i,j)*KXX+SXX0(i,j)*(1-KXX);
            SYY(i,j)=2*ETAP(i,j)*EYY(i,j)*KXX+SYY0(i,j)*(1-KXX);
            
            DSXX(i,j)=SXX(i,j)-SXX0(i,j);
            DSYY(i,j)=SYY(i,j)-SYY0(i,j);
        end
    end
    
    % Compute stress and strain rate invariants and dissipation
    % Process pressure cells
    for i=2:1:Ny
        for j=2:1:Nx
            % EXY term is averaged from four surrounding basic nodes
            EXY2=(EXY(i,j)^2+EXY(i-1,j)^2+EXY(i,j-1)^2+EXY(i-1,j-1)^2)/4;
            
            % FABIAN
            EZX2=(EZX(i,j)^2+EZX(i-1,j)^2+EZX(i,j-1)^2+EZX(i-1,j-1)^2)/4;
            EZY2=(EZY(i,j)^2+EZY(i-1,j)^2+EZY(i,j-1)^2+EZY(i-1,j-1)^2)/4;
            EII(i,j)=(EXX(i,j)^2+EXY2+EZX2+EZY2)^0.5;
            
            % Second strain rate invariant SII
            % SXY term is averaged from four surrounding basic nodes
            SXY2=(SXY(i,j)^2+SXY(i-1,j)^2+SXY(i,j-1)^2+SXY(i-1,j-1)^2)/4;
            
            % FABIAN
            SZX2=(SZX(i,j)^2+SZX(i-1,j)^2+SZX(i,j-1)^2+SZX(i-1,j-1)^2)/4;
            SZY2=(SZY(i,j)^2+SZY(i-1,j)^2+SZY(i,j-1)^2+SZY(i-1,j-1)^2)/4;
            SII(i,j)=(0.5*(SXX(i,j)^2+SYY(i,j)^2)+SXY2+SZX2+SZY2)^0.5;
            
            % Dissipation
            DISXY(i,j)=(SXY(i,j)^2/2/ETA(i,j)+SXY(i-1,j)^2/2/ETA(i-1,j)+SXY(i,j-1)^2/2/ETA(i,j-1)+SXY(i-1,j-1)^2/2/ETA(i-1,j-1))/4;
            DISZX(i,j)=(SZX(i,j)^2/2/ETA(i,j)+SZX(i-1,j)^2/2/ETA(i-1,j)+SZX(i,j-1)^2/2/ETA(i,j-1)+SZX(i-1,j-1)^2/2/ETA(i-1,j-1))/4;
            DISZY(i,j)=(SZY(i,j)^2/2/ETA(i,j)+SZY(i-1,j)^2/2/ETA(i-1,j)+SZY(i,j-1)^2/2/ETA(i,j-1)+SZY(i-1,j-1)^2/2/ETA(i-1,j-1))/4;
            DIS(i,j)=SXX(i,j)^2/2/ETAP(i,j)+SYY(i,j)^2/2/ETAP(i,j)+2*DISXY(i,j)+2*DISZX(i,j)+2*DISZY(i,j);
            
            
            if(i<Ny)
                pt_ave(i,j)  = (pt(i,j)+pt(i+1,j))/2;
                pf_ave(i,j)  = (pf(i,j)+pf(i+1,j))/2;
                PT0_ave      = (PT0(i,j)+PT0(i+1,j))/2;
                PF0_ave      = (PF0(i,j)+PF0(i+1,j))/2;
            else
                pt_ave(i,j)  = pt(i,j);
                pf_ave(i,j)  = pf(i,j);
                PT0_ave      = PT0(i,j);
                PF0_ave      = PF0(i,j);
            end
            
            % Compute elastic and viscous compaction
            VIS_COMP(i,j) = (pt_ave(i,j)-pf_ave(i,j))/(ETAB(i,j)*(1-POR(i,j)));
            % Drained compressibility
            BETTADRAINED=(1/GGGB(i,j)+BETTASOLID)/(1-POR(i,j));
            % Biott-Willis koefficient
            KBW=1-BETTASOLID/BETTADRAINED;
            EL_DECOM(i,j) = BETTADRAINED*(pt_ave(i,j)-PT0_ave-KBW*pf_ave(i,j)+KBW*PF0_ave)/dt;
            
        end
    end
    
    
    % Runge-Kutta velocity, spin array
    vxm=zeros(4,1);
    vym=zeros(4,1);
    spm=zeros(4,1);
    % Move markers by nodal velocity field
    for m=1:-1:marknum
        % Save marker position
        xold=xm(m);
        yold=ym(m);
        for rk=1:1:4
            % vx-velocity interpolation
            % [i,j]--------[i,j+1]
            %   |                |
            %   |    o m         |
            %   |                |
            % [i+1,j]-------[i+1,j+1]
            % Indexes and distances
            j=fix(xm(m)/dx)+1;
            i=fix((ym(m)+dy/2)/dy)+1;
            if(j<1)
                j=1;
            elseif (j>Nx-1)
                j=Nx-1;
            end
            if(i<1)
                i=1;
            elseif (i>Ny)
                i=Ny;
            end
            %Distances
            dxm=(xm(m)-xvx(j))/dx;
            dym=(ym(m)-yvx(i))/dy;
            % Weights
            wtmij=(1-dxm)*(1-dym);
            wtmi1j=(1-dxm)*(dym);
            wtmij1=(dxm)*(1-dym);
            wtmi1j1=(dxm)*(dym);
            % Interpolation
            vxm(rk)=vxs(i,j)*wtmij+vxs(i+1,j)*wtmi1j+vxs(i,j+1)*wtmij1+vxs(i+1,j+1)*wtmi1j1;
            
            % vy-velocity interpolation
            % [i,j]--------[i,j+1]
            %   |                |
            %   |    o m         |
            %   |                |
            % [i+1,j]-------[i+1,j+1]
            % Indexes and distances
            j=fix((xm(m)+dx/2)/dx)+1;
            i=fix(ym(m)/dy)+1;
            if(j<1)
                j=1;
            elseif (j>Nx)
                j=Nx;
            end
            if(i<1)
                i=1;
            elseif (i>Ny-1)
                i=Ny-1;
            end
            %Distances
            dxm=(xm(m)-xvy(j))/dx;
            dym=(ym(m)-yvy(i))/dy;
            % Weights
            wtmij=(1-dxm)*(1-dym);
            wtmi1j=(1-dxm)*(dym);
            wtmij1=(dxm)*(1-dym);
            wtmi1j1=(dxm)*(dym);
            % Interpolation
            vym(rk)=vys(i,j)*wtmij+vys(i+1,j)*wtmi1j+vys(i,j+1)*wtmij1+vys(i+1,j+1)*wtmi1j1;
            
            
            % vz-velocity interpolation (basic node)
            % ESP=1/2(dVy/dx-dVx/dy) interpolation
            % [i,j]--------[i,j+1]
            %   |                |
            %   |    o m         |
            %   |                |
            % [i+1,j]-------[i+1,j+1]
            % Indexes and distances
            j=fix((xm(m))/dx)+1;
            i=fix((ym(m))/dy)+1;
            if(j<1)
                j=1;
            elseif (j>Nx-1)
                j=Nx-1;
            end
            if(i<1)
                i=1;
            elseif (i>Ny-1)
                i=Ny-1;
            end
            
            %Distances
            dxm=(xm(m)-x(j))/dx;
            dym=(ym(m)-y(i))/dy;
            % Weights
            wtmij=(1-dxm)*(1-dym);
            wtmi1j=(1-dxm)*(dym);
            wtmij1=(dxm)*(1-dym);
            wtmi1j1=(dxm)*(dym);
            % Interpolation ESP=1/2(dVy/dx-dVx/dy) for the marker
            spm(rk)=ESP(i,j)*wtmij+ESP(i+1,j)*wtmi1j+ESP(i,j+1)*wtmij1+ESP(i+1,j+1)*wtmi1j1;
            vzm(rk)=vzs(i,j)*wtmij+vzs(i+1,j)*wtmi1j+vzs(i,j+1)*wtmij1+vzs(i+1,j+1)*wtmi1j1;
            
            % Moving between A,B,C,D points
            if(rk<3)
                % Moving A->B and A->C
                xm(m)=xold+vxm(rk)*dt/2;
                ym(m)=yold+vym(rk)*dt/2;
            elseif(rk==3)
                % Moving A->D
                xm(m)=xold+vxm(rk)*dt;
                ym(m)=yold+vym(rk)*dt;
            end
        end
        % Compute effective velocity, rotation rate
        vxeff=(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4))/6;
        vyeff=(vym(1)+2*vym(2)+2*vym(3)+vym(4))/6;
        %speff=(spm(1)+2*spm(2)+2*spm(3)+spm(4))/6;
        speff=spm(1);
        
        % Rotate stress on marker according to its spin
        % Compute amount of rotation from spin rate:
        % Espin=1/2(dvy/dx-dvx/dy) i.e. positive for clockwise rotation
        % (when x axis is directed rightward and y axis is directed downward)
        dspeff=speff*dt;
        % Save old stresses
        msxxold=sxxm(m);
        msyyold=syym(m);
        msxyold=sxym(m);
        % SxyNEW=0.5(Sxx-Syy)*sin(2*Espin*dt)+Sxy*cos(2*Espin*dt)
        sxym(m)=0.5*(msxxold-msyyold)*sin(2*dspeff)+msxyold*cos(2*dspeff);
        % SxxNEW=Sxx*(cos(Espin*dt))^2+Syy*(sin(Espin*dt))^2-Sxy*sin(2*Espin*dt)
        sxxm(m)=msxxold*(cos(dspeff))^2+msyyold*(sin(dspeff))^2-msxyold*sin(2*dspeff);
        % SyyNEW=Sxx*(sin(Espin*dt))^2+Syy*(cos(Espin*dt))^2+Sxy*sin(2*Espin*dt)
        syym(m)=msxxold*(sin(dspeff))^2+msyyold*(cos(dspeff))^2+msxyold*sin(2*dspeff);
        
        % Move markers
        xm(m)=xold+vxeff*dt;
        ym(m)=yold+vyeff*dt;
        
    end
    
    % Update timesum
    timesum=timesum+dt;
    timesumcur(timestep)=timesum;
    timemyr=timesum/(1e+6*365.25*24*3600);
    timeyr=timesum/(365.25*24*3600);
    timehr=timesum/(3600);
    dtcur(timestep)=dt;
    maxvxsmod(timestep)=-1e+30;
    minvxsmod(timestep)=1e+30;
    maxvysmod(timestep)=-1e+30;
    minvysmod(timestep)=1e+30;
    maxvzsmod(timestep)=-1e+30;
    minvzsmod(timestep)=1e+30;
    
    % Vx
    for i=2:1:Ny
        for j=1:1:Nx
            if(RHOX(i,j)>2000)
                maxvxsmod(timestep)=max(maxvxsmod(timestep),vxs(i,j));
                minvxsmod(timestep)=min(minvxsmod(timestep),vxs(i,j));
            end
        end
    end
    % Vy
    for i=1:1:Ny
        for j=2:1:Nx
            if(RHOY(i,j)>2000)
                maxvysmod(timestep)=max(maxvysmod(timestep),vys(i,j));
                minvysmod(timestep)=min(minvysmod(timestep),vys(i,j));
            end
        end
    end
    
    % FABIAN
    % Vz
    for i=1:1:Ny
        for j=1:1:Nx
            if(RHOY(i,j)>2000)
                maxvzsmod(timestep)=max(maxvzsmod(timestep),vzs(i,j));
                minvzsmod(timestep)=min(minvzsmod(timestep),vzs(i,j));
            end
        end
    end
    
    % Update VX0
    VX0=vxs;
    VY0=vys;
    
    % FABIAN
    VZ0=vzs;
    for i=1:1:Ny1
        for j=1:1:Nx1
            if(POR(i,j)>0)
                VXF0(i,j)=vxs(i,j)+vxD(i,j)/POR(i,j);
                VYF0(i,j)=vys(i,j)+vyD(i,j)/POR(i,j);
            end
        end
    end
    
    % Update SXX0
    SXX0=SXX;
    % Update SYY0
    SYY0=SYY;
    % Update SXY0
    SXY0=SXY;
    
    % FABIAN
    % Update SZX0
    SZX0=SZX;
    % Update SZY0
    SZY0=SZY;
    % Update PTF0
    PTF0=(pt-pf);
    PT0=pt;
    PF0=pf;
    
    
    fprintf('==================================== \n');
    fprintf('total time:        %.12E sec \n',timesum);
    fprintf('time step:         %.12E sec \n',dt);
    fprintf('Vslip max:         %.12E m/s \n',Vmax);
    fprintf('iter-iterations:   %d \n',iterstep);
    fprintf('global-iterations: %d \n',ynlast);
        
    if seismic_cycles_data
        
        if (timesum==dt)
            data_save = x;
            fileID    = fopen('x_fault.txt','a');
            TP_write  = fprintf (fileID,'%.3E    \n',data_save);
            fclose(fileID);
        end
        
        if timesum_plus < timesum
            
            % ========== save slip rate
            data_save = [timesum dt VSLIPB(line_fault,:)];
            fid = fopen('EVO_Vslip.txt','a');
            fprintf(fid,'%.6E ',data_save);
            fprintf(fid,'\n');
            fclose(fid);
            clear data_save
            
            % ========== save viscosity
            data_save = [timesum dt ETA(line_fault,:)];
            fid = fopen('EVO_viscosity.txt','a');
            fprintf(fid,'%.6E ',data_save);
            fprintf(fid,'\n');
            fclose(fid);
            clear data_save
            
            % ========== save fluid pressure
            data_save = [timesum dt pf(line_fault,:)];
            fid = fopen('EVO_press_flu.txt','a');
            fprintf(fid,'%.6E ',data_save);
            fprintf(fid,'\n');
            fclose(fid);
            clear data_save
            
            % ========== save effective pressure
            P_diff = pt-pf;
            data_save = [timesum dt P_diff(line_fault,:)];
            fid = fopen('EVO_press_eff.txt','a');
            fprintf(fid,'%.6E ',data_save);
            fprintf(fid,'\n');
            fclose(fid);
            clear data_save
            
            % ========== save SigmaY
            data_save = [timesum dt SigmaY(line_fault,:)];
            fid = fopen('EVO_SigmaY.txt','a');
            fprintf(fid,'%.9E ',data_save);
            fprintf(fid,'\n');
            fclose(fid);
            clear data_save
            
            % ========== save SII
            data_save = [timesum dt SII_fault(line_fault,:)];
            fid = fopen('EVO_Sii.txt','a');
            fprintf(fid,'%.9E ',data_save);
            fprintf(fid,'\n');
            fclose(fid);
            clear data_save
            
            % ========== save Theta
            data_save = [timesum dt OM(line_fault,:)];
            fid = fopen('EVO_Theta.txt','a');
            fprintf(fid,'%.9E ',data_save);
            fprintf(fid,'\n');
            fclose(fid);
            clear data_save
            
            % ========== save viscous compaction
            data_save = [timesum dt VIS_COMP(line_fault,:)];
            fid = fopen('EVO_Visc_comp.txt','a');
            fprintf(fid,'%.9E ',data_save);
            fprintf(fid,'\n');
            fclose(fid);
            clear data_save
            
            % ========== save elastic compaction
            data_save = [timesum dt EL_DECOM(line_fault,:)];
            fid = fopen('EVO_Elast_comp.txt','a');
            fprintf(fid,'%.9E ',data_save);
            fprintf(fid,'\n');
            fclose(fid);
            clear data_save
            
            % ========== save time, dt, vmax
            Vmax = max(max(VSLIPB));
            fileID    = fopen('EVO_data.txt','a');
            TP_write  = fprintf (fileID,'%.12E  %.12E  %.12E  %d   %d \n',timesum,dt,Vmax,ynlast,iterstep);
            fclose(fileID);
            
            timesum_plus = timesum;
        end
    end
    
    if(fix(timestep/savematstep)*savematstep==timestep)
        namemat    =  [nname,num2str(timestep)];
        save(namemat);
        fdata=fopen('file.txt','wt');
        fprintf(fdata,'%d',timestep);
        fclose(fdata);
    end
            
end
