% ======================================== 
% H - MEC: Hydro - Mechanical Earthquake Cycle
% Computational Earthquake Physics
% ETH Zurich,  2022
%
% Dal Zilio,  L.,  Hegyi,  B.,  Behr,  W. M.,  & Gerya,  T. (2022)
% Hydro - mechanical earthquake cycles in a
% poro - visco - elasto - plastic fluid - bearing fault.
% DOI: http://arxiv.org / abs / 2201.11786
% ======================================== 
% Solving of compressible Stokes + continuity equations
% with variable viscosity and Darsi + continuity equations
% Total X - Stokes: dSIGMAxxt' / dx + dSIGMAxyt' / dy - dPt / dx = -RHOt * gx
% Total Y - Stokes: dSIGMAyxt' / dx + dSIGMAyyt' / dy - dPt / dy = -RHOt * gy
% Solid Continuity: dVxs / dx + dVys / dy + (Pt - Pf) / ETAbulk = 0
% Fluid X - Darsi:  - ETAfluid / K * VxD - dPf / dx = -RHOf * gx
% Fluid Y - Darsi:  - ETAfluid / K * VyD - dPf / dy = -RHOf * gy
% Fluid Continuity: dVxD / dx + dVyD / dy - (Pt - Pf) / ETAbulk = 0
% + staggered grid
% + P - v formulation
% + advecting material fileds by markers
% ======================================== 

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% to do:
%
% 1. Change variable names to useful / understandable names
%
% 2. Move all variable declarations outside of loops if possible
%
% 3. Do calculations outside of loops if possible
%   a) Multiplications of every element can be done at once outside of loop
%
% 4. Remove unnessesary calculations if there are any
%   a) This includes copying variables not needed
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% /////////////////////////////////////////////////////////////////////////////////////// 
% inputs
% /////////////////////////////////////////////////////////////////////////////////////// 

warning off
% Load mat file
fdata = fopen('file.txt', 'rt');
timestep = fscanf(fdata, '%d', 1);
fclose(fdata);

num_timesteps = 10; % set timestep limit

if(timestep  >  0)
    namemat = ['h_mec_',  num2str(timestep)];
    load(namemat);
    timestep = timestep + 1;
else
    % ======================================== 
    % Start new run
    % (1) Clearing
    clear all;
    
    % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    % start measuring time
    t_start = tic;
    t_init = tic;
    % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    % /////////////////////////////////////////////////////////////////////////////////////// 
    % global variable declarations
    % /////////////////////////////////////////////////////////////////////////////////////// 
    
    % (2) Define Numerical model
    % Eulerian basic grid
    xsize = 5000.0;         % size in horizontal direction,  m
    ysize = 1250.0;         % size in vertical direction,  m
    Nx = 21;                % number of grid steps in horizontal directions
    Ny = 6;                 % number of grid steps in vertical direction
    % Eulerian Staggered Grid
    Nx1 = Nx + 1;              % Number of horizontal lines for staggered grid
    Ny1 = Ny + 1;              % Number of vertical lines for staggered grid
    
    % ======================================== 
    % Output files
    seismic_cycles_data = 1;
    timesum_plus        = 0;
    line_fault          = (Ny - 1)  /  2 + 1;
    % ======================================== 
    % Coordinates
    dx = xsize / (Nx - 1); % grid step in horizontal direction,  m
    dy = ysize / (Ny - 1); % grid step in vertical direction,  m
    xbeg = 0;
    xend = xsize;
    ybeg = 0;
    yend = ysize;
    x = xbeg:dx:xend; % horizontal coordinates of basic grid points
    y = ybeg:dy:yend; % vertical coordinates of basic grid points
    xvx = xbeg:dx:xend + dx; % Horizontal coordinates of Vx - nodes
    yvx = ybeg - dy / 2:dy:yend + dy / 2; % Vertical coordinates of Vx - nodes
    xvy = xbeg - dx / 2:dx:xend + dx / 2; % Horizontal coordinates of Vy - nodes
    yvy = ybeg:dy:yend + dy; % Vertical coordinates of Vy - nodes
    xp = xbeg - dx / 2:dx:xend + dx / 2; % Horizontal coordinates of P - nodes
    yp = ybeg - dy / 2:dy:yend + dy / 2; % Vertical coordinates of Vx - nodes
    
    %          Block     Fault
    arsfm =   [0.0300    0.0090 ]; % a - parameter of RSF
    brsfm =   [0.0010    0.0150 ]; % b - parameter of RSF
    lrsfm =   [0.0140    0.0140 ]; % L - parameter of RSF (characteristic slip distance)
    omm   =   [10       -5      ]; % State
    V0    =   1e-9;                % Reference slip velocity of RSF,  m / s
    alpha =   29;                  % Scaling factor for shear viscosity of porous matrix, 
    % ======================================== 
    POR0 = 0.01; % Standard porosity
    % ======================================== 
    % Brittle / plastic rheology
    cohes = 0.0e+6;     % Cohesion,  Pa
    friction = 0.6;     % Internal friction coefficient confined
    dilatation = 0.00;  % Dilatation coefficient confined
    tensile = 1.0;      % Internal friction coefficient tensile
    shearmod = 3.0e+10; % Shear modulus
    BETTAFLUID = 1e-8;  % 4.0e-10; % Compressibility of fluid,  1 / Pa
    BETTASOLID = 2e-11; % 2.5e-11; % Compressibility of solid,  1 / Pa
    PORDIL = 0.15;      % Upper limit of porosity for dilatation
    ESLIP0 = 1e+04;     % ETHAslip LN model,  Pa s
    VSLIP0 = 0.5e-9;    % Characteristic velocity,  m / s
    faultwidth = dx;    % Characteristic fault width,  m
    % ======================================== 
    % Constants
    gx = 0;             % Horizontal gravity,  m / s^2
    gy = 0;             % Vertical gravity,  m / s^2
    PCONF = 10.0e6;     % Confining pressure
    PTFDIFF = 30e+6;    % Total - Fluid Pressure difference in the top row,  Pa
    SIGMALOAD = 0;      % SIGMAyy at the top,  Pa
    ETATOP = 1e-3;      % Viscosity of the sticky water,  Pa * s
    % ======================================== 
    % Limits
    pormin = 1e-4; % Min porosity limit
    pormax = 1 - pormin; % Max porosity limit
    etamin = 1e-3; % Lower shear viscosity cutoff
    etamax = 1e+50; % Upper shear viscosity cutoff
    etabmin = 1e-3; % Lower bulk viscosity cutoff
    etabmax = 1e+25; % Upper bulk viscosity cutoff
    gmin = 1e+08; % Lower shear modulus cutoff
    gmax = 1e+11; % Upper shear modulus cutoff
    kkkmin = 1e-22; % Lower Darsi viscosity cutoff
    kkkmax = 1e-12; % Upper Darsi viscosity cutoff
    % ======================================== 
    % Limits
    dsubgrids = 0; % Subgrid diffusion for stresses
    dsubgridv = 0; % Subgrid diffusion for velocity
    stpmax = 2e-4;% / dy * faultwidth; % Max gridstep fraction for marker displacement in the channel
    stpmax1 = 6e-5;% / dy * faultwidth; % Max gridstep fraction for marker displacement
    peffmax = PTFDIFF; % Pressure to normalise timestep
    % ======================================== 
    
    % Basic nodes
    OM0 = ones(Ny, Nx);    % Old state parameter
    OM = ones(Ny, Nx);     % State parameter
    ARSF = ones(Ny, Nx); % a - parameter of RSF
    BRSF = ones(Ny, Nx); % b - parameter of RSF
    LRSF = ones(Ny, Nx); % L - parameter of RSF
    
    % Unknown parameters
    pt = zeros(Ny1, Nx1);  % Total pressure
    vxs = zeros(Ny1, Nx1); % Solid vx - velocity
    vys = zeros(Ny1, Nx1); % Solid vy - velocity
    pf = zeros(Ny1, Nx1);  % Fluid pressure
    vxD = zeros(Ny1, Nx1); % Darsi vx - velocity
    vyD = zeros(Ny1, Nx1); % Darsi vy - velocity
    
    % Nodal arrays
    % Basic nodes
    RHO = zeros(Ny, Nx);
    ETA = zeros(Ny, Nx);
    ETA0 = zeros(Ny, Nx);
    IETAPLB = zeros(Ny, Nx);
    SXY = zeros(Ny, Nx);
    SXY0 = zeros(Ny, Nx);
    YNY0 = zeros(Ny, Nx);
    KKK = zeros(Ny, Nx);
    GGG = zeros(Ny, Nx);
    COHC = zeros(Ny, Nx);
    COHT = zeros(Ny, Nx);
    FRIC = zeros(Ny, Nx);
    FRIT = zeros(Ny, Nx);
    DILC = zeros(Ny, Nx);
    TTT = zeros(Ny, Nx);
    EIIB = zeros(Ny, Nx);
    STRPLB = zeros(Ny, Nx);
    
    lower_block = ysize / 2 + (dy);
    upper_block = ysize / 2 - (dy);
    
    BRSF(:, :) = brsfm(1);
    ARSF(:, :) = arsfm(1);
    LRSF(:, :) = lrsfm(1);
    
    
    % ======================== 
    % Where to apply the transition on the left(1) and right(2)
    TS_1 = 6e3;
    TS_2 = 8e3;
    
    TS_3 = 34e3;
    TS_4 = 37e3;
    
    % Define Fault
    for j = 1:1:Nx
        for i = 1:1:Ny
            
            OM0(i, j) = omm(1);
            
            if (y(i) > upper_block && y(i) < lower_block)
                OM0(i, j) = omm(2);
            end
            
            if(x(j) < TS_1 && y(i) > upper_block && y(i) < lower_block)
                BRSF(i, j) = brsfm(1);
                ARSF(i, j) = arsfm(1);
            end
            if(x(j) >= TS_1 && x(j) < TS_2 && y(i) > upper_block && y(i) < lower_block)
                BRSF(i, j) = brsfm(1) - (brsfm(1) - brsfm(2)) * ((x(j) - TS_1) / (TS_2 - TS_1));
                ARSF(i, j) = arsfm(1) - (arsfm(1) - arsfm(2)) * ((x(j) - TS_1) / (TS_2 - TS_1));
            end
            if(x(j) >= TS_2 && x(j) <= TS_3 && y(i) > upper_block && y(i) < lower_block)
                BRSF(i, j) = brsfm(2);
                ARSF(i, j) = arsfm(2);
            end
            if(x(j) > TS_3 && x(j) <= TS_4 && y(i) > upper_block && y(i) < lower_block)
                BRSF(i, j) = brsfm(2) - (brsfm(2) - brsfm(1)) * ((x(j) - TS_3) / (TS_4 - TS_3));
                ARSF(i, j) = arsfm(2) - (arsfm(2) - arsfm(1)) * ((x(j) - TS_3) / (TS_4 - TS_3));
            end
            if(x(j) > TS_4 && y(i) > upper_block && y(i) < lower_block)
                BRSF(i, j) = brsfm(1);
                ARSF(i, j) = arsfm(1);
            end
        end
    end
    
    OM = OM0; % ??? not sure if necessary
    
    % Pressure nodes
    ETAB = zeros(Ny1, Nx1);
    ETAB0 = zeros(Ny1, Nx1);
    ETAP = zeros(Ny1, Nx1);
    ETAP0 = zeros(Ny1, Nx1);
    POR = zeros(Ny1, Nx1);
    GGGP = zeros(Ny1, Nx1);
    GGGB = zeros(Ny1, Nx1);
    PTF0 = zeros(Ny1, Nx1);
    PT0 = (PCONF + PTFDIFF) * ones(Ny1, Nx1);
    PF0 = PCONF * ones(Ny1, Nx1);
    pt_ave = zeros(Ny1, Nx1);
    pf_ave = zeros(Ny1, Nx1);
    SXX = zeros(Ny1, Nx1);
    SXX0 = zeros(Ny1, Nx1);
    SYY = zeros(Ny1, Nx1);
    SYY0 = zeros(Ny1, Nx1);
    DILP = zeros(Ny1, Nx1);
    % Vx nodes
    RHOX = zeros(Ny1, Nx1);
    RHOFX = zeros(Ny1, Nx1);
    ETADX = zeros(Ny1, Nx1);
    PORX = zeros(Ny1, Nx1);
    VX0 = zeros(Ny1, Nx1);
    VXF0 = zeros(Ny1, Nx1);
    % Vy nodes
    RHOY = zeros(Ny1, Nx1);
    RHOFY = zeros(Ny1, Nx1);
    ETADY = zeros(Ny1, Nx1);
    PORY = zeros(Ny1, Nx1);
    VY0 = zeros(Ny1, Nx1);
    VYF0 = zeros(Ny1, Nx1);
    
    VSLIPB = zeros(Ny, Nx);
    
    
    % Lagrangian solid markers
    Nxm = (Nx - 1) * 4; % Marker resolution in x - dection
    Nym = (Ny - 1) * 4; % Marker resolution in y direction
    dxms = xsize / Nxm; % Standard marker horizontal step
    dyms = ysize / Nym; % Standard marker vertical step
    marknum = Nxm * Nym; % Total number of markers
    rhom = zeros(marknum, 1); % Density of solid
    etasm = zeros(marknum, 1); % Standard shear viscosity of bulk
    etam = zeros(marknum, 1); % Shear viscosity of bulk
    etabm = zeros(marknum, 1); % Bulk viscosity of bulk
    cohescm = zeros(marknum, 1); % Cohesion for confined fracture of solid
    frictcm = zeros(marknum, 1); % friction for confined fracture of solid
    dilatcm = zeros(marknum, 1); % dilatation for confined fracture of solid
    cohestm = zeros(marknum, 1); % Cohesion for tensile fracture of solid
    fricttm = zeros(marknum, 1); % friction for tensile fracture of solid
    porm = zeros(marknum, 1); % Porosity of solid
    kkkm = zeros(marknum, 1); % Standard permeability of solid
    rhofm = zeros(marknum, 1); % Density of fluid
    etafm = zeros(marknum, 1); % Viscosity of fluid
    tm = zeros(marknum, 1); % Marker rock type
    xm = zeros(marknum, 1); % Horizontal coordinates of solid markers
    ym = zeros(marknum, 1); % Vertical coordinates of solid markers
    sxxm = zeros(marknum, 1); % Marker SIGMAxx',  Pa
    syym = zeros(marknum, 1); % Marker SIGMAyy',  Pa
    sxym = zeros(marknum, 1); % Marker SIGMAxy',  Pa
    gsm = zeros(marknum, 1); % Standard shear modulus of bulk,  Pa
    gm = zeros(marknum, 1); % Shear modulus of bulk,  Pa
    vx0m = zeros(marknum, 1); % Marker horizontal velocity
    vy0m = zeros(marknum, 1); % Marker vertical velocity
    ptfm = zeros(marknum, 1); % Pt - Pf,  Pa
    amursfm = zeros(marknum, 1); % RSF a / mu parameter
    
    
    m = 1;
    for jm = 1:1:Nxm
        for im = 1:1:Nym
            % Define randomized regular coordinates
            xm(m) = xbeg + dxms  /  2 + (jm - 1) * dxms + (rand - 0.5) * dxms;
            ym(m) = ybeg + dyms  /  2 + (im - 1) * dyms + (rand - 0.5) * dyms;
            %Matrix
            tm(m) = 1;
            rhom(m) = 3000;
            etasm(m) = 10^22;
            gsm(m) = shearmod;
            cohescm(m) = cohes;
            cohestm(m) = cohes;
            frictcm(m) = friction;
            dilatcm(m) = dilatation;
            fricttm(m) = tensile;
            porm(m) = 0.01 * (1 + 0.0 * (rand - 0.5));
            kkkm(m) = 5e-16;% * (dy / faultwidth)^2;
            rhofm(m) = 1000;
            etafm(m) = 1e-3;
            etam(m) = etasm(m) * exp( - alpha * porm(m));
            etabm(m) = etam(m) / porm(m);% / (1 - porm(m));
            gm(m) = gsm(m);% * (1 - porm(m));
            
            % Air,  wedge,  slab
            if(ym(m) < upper_block || ym(m) > lower_block)
                tm(m) = -1;
                etam(m) = 1e+23;
                etasm(m) = 1e+23;
                rhom(m) = 3000;
                kkkm(m) = 5e-16;% * (dy / faultwidth)^2;
                cohescm(m) = cohes * 1e+3;
                cohestm(m) = cohes * 1e+3;
                frictcm(m) = friction;
                dilatcm(m) = dilatation;
                fricttm(m) = tensile;
            end
            
            % Update marker index
            m = m + 1;
        end
    end
    % Average porosity
    poravr0 = mean(porm);
    
    
    % (3) Defining global matrixes
    % according to the global number of unknowns
    N = Nx1 * Ny1 * 6; % Global number of unknowns
    L = sparse(N, N); % Matrix of coefficients in the left part
    R = zeros(N, 1); % Vector of the right parts of equations
    
    
    % 4) Boundary conditions
    bcleft = 1;
    bcright = 1;
    bcupper = -2e-9;
    bclower = 2e-9;
    bcvyflower = 0;% - 1e-12;
    % bcvxfleft = bcvyflower * xsize / ysize;
    
    
    % Entering time cycle
    timesum = 0;
    dtelastic0 = 5e+8; % elastic timestep
    dt = dtelastic0;
    dtmin = 1e-4;
    
    ascale = 1e+0;
    dtreset = 1e-20; % Minimum timestep for velocity reset
    
    % /////////////////////////////////////////////////////////////////////////////////////// 
    % move to global variables
    % /////////////////////////////////////////////////////////////////////////////////////// 
    
    % Timesteps between visualization frames
    savematstep = 300;  % .mat storage periodicity
    nname  = 'h_mec_';  % mat filename
    niterglobal = 10000; % Max number of global iterations
    ynlastmax = 499; % Max number of iterations without changing timestep
    dtstep = 5; % Min number of plastic iterations before changing dt
    vratiomax = 0.001; % Max change of velocity between iterations
    dtkoef = 1.1; % Koefficient for decrease of previous timestep
    dtkoefup = 1.2; % Koefficient for increase of previous timestep
    dtkoefv = 1.001; % koef for velocity - based timestep reduction
    errmin = 1e+2; % min LSQ err for stoping iterations
    etawt = 0.0;
    syieldmin = 1e-3;
    tyield = 1;
    timestep = 1;
    maxvxsmod = 0;
    minvxsmod = 0;
    maxvysmod = 0;
    minvysmod = 0;
    yndtdecrease = 1;
    plstrain = 0; % Plastic strain for dilatancy
end

aa1 = 20.93;
aa2 = 35.28;
aa3 = 2.34;
bb1 = 0.99;
bb2 = 44.39;
bb3 = 0.73;
cc1 = 0.37;
cc2 = 3.54;
cc3 = 0.47;
% Plastic Strain:
gammapij = plstrain * 200;
gammapi1j = plstrain * 200;
gammapij1 = plstrain * 200;
gammapi1j1 = plstrain * 200;

% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
% added these variables for time measurement purpose
total_time_markers = 0;
total_nodes = 0;
total_time_global = 0;
total_plastic_strain = 0;
total_viscoplastic = 0;
total_global_matrices = 0;
total_solve = 0;
total_reload_solution = 0;
total_internal_nodes = 0;
total_pressure_cell = 0;
total_pressure_cell2 = 0;
total_plasticity = 0;
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
% stop measuring initialisation time
initialisation_time = toc(t_init);
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

% /////////////////////////////////////////////////////////////////////////////////////// 
% actual computations start here
% /////////////////////////////////////////////////////////////////////////////////////// 

% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
% start measuring iteration time
t_iteration = tic;
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

num_timesteps = 10;

for timestep = timestep:1:num_timesteps
    
    % /////////////////////////////////////////////////////////////////////////////////////// 
    % declare variables outside of loop and set zero at the start of each step
    % /////////////////////////////////////////////////////////////////////////////////////// 
    
    % Interpolate ETA,  RHO to nodal points
    % Basic nodes
    RHOSUM = zeros(Ny, Nx);
    ETASUM = zeros(Ny, Nx);
    KKKSUM = zeros(Ny, Nx);
    TTTSUM = zeros(Ny, Nx);
    SXYSUM = zeros(Ny, Nx);
    GGGSUM = zeros(Ny, Nx);
    ETA0SUM = zeros(Ny, Nx);
    COHCSUM = zeros(Ny, Nx);
    FRICSUM = zeros(Ny, Nx);
    DILCSUM = zeros(Ny, Nx);
    COHTSUM = zeros(Ny, Nx);
    FRITSUM = zeros(Ny, Nx);
    WTSUM = zeros(Ny, Nx);
    
    %LDZ
    OM0SUM = zeros(Ny, Nx); % Old state parameter
    OMSUM = zeros(Ny, Nx); % State parameter
    ARSFSUM = zeros(Ny, Nx); % a - parameter of RSF
    BRSFSUM = zeros(Ny, Nx); % b - parameter of RSF
    LRSFSUM = zeros(Ny, Nx); % L - parameter of RSF
    
    % Pressure nodes
    ETAPSUM = zeros(Ny1, Nx1);
    ETAP0SUM = zeros(Ny1, Nx1);
    ETAB0SUM = zeros(Ny1, Nx1);
    PORSUM = zeros(Ny1, Nx1);
    SXXSUM = zeros(Ny1, Nx1);
    SYYSUM = zeros(Ny1, Nx1);
    GGGPSUM = zeros(Ny1, Nx1);
    WTPSUM = zeros(Ny1, Nx1);
    % Vx nodes
    RHOXSUM = zeros(Ny1, Nx1);
    RHOFXSUM = zeros(Ny1, Nx1);
    ETADXSUM = zeros(Ny1, Nx1);
    PORXSUM = zeros(Ny1, Nx1);
    VX0SUM = zeros(Ny1, Nx1);
    WTXSUM = zeros(Ny1, Nx1);
    % Vy nodes
    RHOYSUM = zeros(Ny1, Nx1);
    RHOFYSUM = zeros(Ny1, Nx1);
    ETADYSUM = zeros(Ny1, Nx1);
    PORYSUM = zeros(Ny1, Nx1);
    VY0SUM = zeros(Ny1, Nx1);
    WTYSUM = zeros(Ny1, Nx1);
    
    % Cycle on markers
    
    % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    % start measuring marker iteration time
    t_marker = tic;
    % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    % /////////////////////////////////////////////////////////////////////////////////////// 
    % second layer loop (timestep -> markers) possibly very costly
    % /////////////////////////////////////////////////////////////////////////////////////// 
    
    for m = 1:1:marknum
        
        % Marker properties
        kkkmm = kkkm(m) * (porm(m) / POR0)^3;
        % Checking permeability limits
        if(kkkmm < kkkmin)
            kkkmm = kkkmin;
        end
        if(kkkmm > kkkmax)
            kkkmm = kkkmax;
        end
        
        % Cohesion,  friction of porous matrix
        
        % Viscosity of porous matrix
        if(tm(m) ~= 0)
            cohescmm = cohescm(m) * (1 - porm(m)); %  * exp( - alpha * porm(m));
            cohestmm = cohestm(m) * (1 - porm(m)); %  * exp( - alpha * porm(m));
            frictcmm = frictcm(m);
            dilatcmm = dilatcm(m);
            fricttmm = fricttm(m);
            etasmm0 = etasm(m);
            etamm0 = etasm(m) * exp( - alpha * porm(m));
            etamm = etam(m);
            % total density
            rhomm = rhom(m) * (1 - porm(m)) + rhofm(m) * porm(m);
        else
            cohescmm = cohescm(m);
            cohestmm = cohestm(m);
            frictcmm = frictcm(m);
            dilatcmm = 0;
            fricttmm = fricttm(m);
            etasmm0 = etasm(m);
            etamm0 = etasm(m);
            etamm = etam(m);
            % total density
            rhomm = rhom(m);
        end
        % Matrix viscosity
        if(etamm0 < etamin)
            etamm0 = etamin;
        end
        if(etamm0 > etamax)
            etamm0 = etamax;
        end
        % Effective viscosity
        if(etamm < etamin)
            etamm = etamin;
        end
        if(etamm > etamax)
            etamm = etamax;
        end
        
        % Darsi "viscosity"
        etadm = etafm(m) / kkkmm;
        
        % Interpolate to basic nodes
        % [i, j] -------- [i, j + 1]
        %   |                |
        %   |    o m         |
        %   |                |
        % [i + 1, j] ------- [i + 1, j + 1]
        % Indexes and distances
        j = fix(xm(m) / dx) + 1;
        i = fix(ym(m) / dy) + 1;
        if(j < 1)
            j = 1;
        elseif(j > Nx - 1)
            j = Nx - 1;
        end
        if(i < 1)
            i = 1;
        elseif(i > Ny - 1)
            i = Ny - 1;
        end
        dxm = (xm(m) - x(j)) / dx;
        dym = (ym(m) - y(i)) / dy;
        
        % /////////////////////////////////////////////////////////////////////////////////////// 
        % many very similar computations. Maybe write a function to do it instead
        % /////////////////////////////////////////////////////////////////////////////////////// 
        
        wtmij = (1 - dxm) * (1 - dym);
        ETASUM(i, j) = ETASUM(i, j) + etamm * wtmij;
        RHOSUM(i, j) = RHOSUM(i, j) + rhomm * wtmij;
        KKKSUM(i, j) = KKKSUM(i, j) + kkkmm * wtmij;
        TTTSUM(i, j) = TTTSUM(i, j) + tm(m) * wtmij;
        SXYSUM(i, j) = SXYSUM(i, j) + sxym(m) * wtmij;
        GGGSUM(i, j) = GGGSUM(i, j) + 1 / gm(m) * wtmij;
        ETA0SUM(i, j) = ETA0SUM(i, j) + etamm0 * wtmij;
        COHTSUM(i, j) = COHTSUM(i, j) + cohestmm * wtmij;
        FRITSUM(i, j) = FRITSUM(i, j) + fricttmm * wtmij;
        COHCSUM(i, j) = COHCSUM(i, j) + cohescmm * wtmij;
        FRICSUM(i, j) = FRICSUM(i, j) + frictcmm * wtmij;
        DILCSUM(i, j) = DILCSUM(i, j) + dilatcmm * wtmij;
        
        WTSUM(i, j) = WTSUM(i, j) + wtmij;
        
        wtmi1j = (1 - dxm) * (dym);
        ETASUM(i + 1, j) = ETASUM(i + 1, j) + etamm * wtmi1j;
        RHOSUM(i + 1, j) = RHOSUM(i + 1, j) + rhomm * wtmi1j;
        KKKSUM(i + 1, j) = KKKSUM(i + 1, j) + kkkmm * wtmi1j;
        TTTSUM(i + 1, j) = TTTSUM(i + 1, j) + tm(m) * wtmi1j;
        SXYSUM(i + 1, j) = SXYSUM(i + 1, j) + sxym(m) * wtmi1j;
        GGGSUM(i + 1, j) = GGGSUM(i + 1, j) + 1 / gm(m) * wtmi1j;
        ETA0SUM(i + 1, j) = ETA0SUM(i + 1, j) + etamm0 * wtmi1j;
        COHTSUM(i + 1, j) = COHTSUM(i + 1, j) + cohestmm * wtmi1j;
        FRITSUM(i + 1, j) = FRITSUM(i + 1, j) + fricttmm * wtmi1j;
        COHCSUM(i + 1, j) = COHCSUM(i + 1, j) + cohescmm * wtmi1j;
        FRICSUM(i + 1, j) = FRICSUM(i + 1, j) + frictcmm * wtmi1j;
        DILCSUM(i + 1, j) = DILCSUM(i + 1, j) + dilatcmm * wtmi1j;
        
        WTSUM(i + 1, j) = WTSUM(i + 1, j) + wtmi1j;
        
        wtmij1 = (dxm) * (1 - dym);
        ETASUM(i, j + 1) = ETASUM(i, j + 1) + etamm * wtmij1;
        RHOSUM(i, j + 1) = RHOSUM(i, j + 1) + rhomm * wtmij1;
        KKKSUM(i, j + 1) = KKKSUM(i, j + 1) + kkkmm * wtmij1;
        TTTSUM(i, j + 1) = TTTSUM(i, j + 1) + tm(m) * wtmij1;
        SXYSUM(i, j + 1) = SXYSUM(i, j + 1) + sxym(m) * wtmij1;
        GGGSUM(i, j + 1) = GGGSUM(i, j + 1) + 1 / gm(m) * wtmij1;
        ETA0SUM(i, j + 1) = ETA0SUM(i, j + 1) + etamm0 * wtmij1;
        COHTSUM(i, j + 1) = COHTSUM(i, j + 1) + cohestmm * wtmij1;
        FRITSUM(i, j + 1) = FRITSUM(i, j + 1) + fricttmm * wtmij1;
        COHCSUM(i, j + 1) = COHCSUM(i, j + 1) + cohescmm * wtmij1;
        FRICSUM(i, j + 1) = FRICSUM(i, j + 1) + frictcmm * wtmij1;
        DILCSUM(i, j + 1) = DILCSUM(i, j + 1) + dilatcmm * wtmij1;
        
        WTSUM(i, j + 1) = WTSUM(i, j + 1) + wtmij1;
        
        wtmi1j1 = (dxm) * (dym);
        ETASUM(i + 1, j + 1) = ETASUM(i + 1, j + 1) + etamm * wtmi1j1;
        RHOSUM(i + 1, j + 1) = RHOSUM(i + 1, j + 1) + rhomm * wtmi1j1;
        KKKSUM(i + 1, j + 1) = KKKSUM(i + 1, j + 1) + kkkmm * wtmi1j1;
        TTTSUM(i + 1, j + 1) = TTTSUM(i + 1, j + 1) + tm(m) * wtmi1j1;
        SXYSUM(i + 1, j + 1) = SXYSUM(i + 1, j + 1) + sxym(m) * wtmi1j1;
        GGGSUM(i + 1, j + 1) = GGGSUM(i + 1, j + 1) + 1 / gm(m) * wtmi1j1;
        ETA0SUM(i + 1, j + 1) = ETA0SUM(i + 1, j + 1) + etamm0 * wtmi1j1;
        COHTSUM(i + 1, j + 1) = COHTSUM(i + 1, j + 1) + cohestmm * wtmi1j1;
        FRITSUM(i + 1, j + 1) = FRITSUM(i + 1, j + 1) + fricttmm * wtmi1j1;
        COHCSUM(i + 1, j + 1) = COHCSUM(i + 1, j + 1) + cohescmm * wtmi1j1;
        FRICSUM(i + 1, j + 1) = FRICSUM(i + 1, j + 1) + frictcmm * wtmi1j1;
        DILCSUM(i + 1, j + 1) = DILCSUM(i + 1, j + 1) + dilatcmm * wtmi1j1;
        
        WTSUM(i + 1, j + 1) = WTSUM(i + 1, j + 1) + wtmi1j1;
        %     end
        
        % Interpolate to pressure nodes
        % [i, j] -------- [i, j + 1]
        %   |                |
        %   |    o m         |
        %   |                |
        % [i + 1, j] ------- [i + 1, j + 1]
        % Indexes and distances
        j = fix((xm(m) + dx / 2) / dx) + 1;
        i = fix((ym(m) + dy / 2) / dy) + 1;
        if(j < 1)
            j = 1;
        elseif(j > Nx)
            j = Nx;
        end
        if(i < 1)
            i = 1;
        elseif(i > Ny)
            i = Ny;
        end
        
        dxm = (xm(m) - xp(j)) / dx;
        dym = (ym(m) - yp(i)) / dy;
        
        wtmij = (1 - dxm) * (1 - dym);
        ETAPSUM(i, j) = ETAPSUM(i, j) + etamm * wtmij;
        PORSUM(i, j) = PORSUM(i, j) + porm(m) * wtmij;
        SXXSUM(i, j) = SXXSUM(i, j) + sxxm(m) * wtmij;
        SYYSUM(i, j) = SYYSUM(i, j) + syym(m) * wtmij;
        GGGPSUM(i, j) = GGGPSUM(i, j) + 1 / gm(m) * wtmij;
        ETAP0SUM(i, j) = ETAP0SUM(i, j) + etamm0 * wtmij;
        ETAB0SUM(i, j) = ETAB0SUM(i, j) + etasmm0 * wtmij;
        WTPSUM(i, j) = WTPSUM(i, j) + wtmij;
        
        wtmi1j = (1 - dxm) * (dym);
        ETAPSUM(i + 1, j) = ETAPSUM(i + 1, j) + etamm * wtmi1j;
        PORSUM(i + 1, j) = PORSUM(i + 1, j) + porm(m) * wtmi1j;
        SXXSUM(i + 1, j) = SXXSUM(i + 1, j) + sxxm(m) * wtmi1j;
        SYYSUM(i + 1, j) = SYYSUM(i + 1, j) + syym(m) * wtmi1j;
        GGGPSUM(i + 1, j) = GGGPSUM(i + 1, j) + 1 / gm(m) * wtmi1j;
        ETAP0SUM(i + 1, j) = ETAP0SUM(i + 1, j) + etamm0 * wtmi1j;
        ETAB0SUM(i + 1, j) = ETAB0SUM(i + 1, j) + etasmm0 * wtmi1j;
        WTPSUM(i + 1, j) = WTPSUM(i + 1, j) + wtmi1j;
        
        wtmij1 = (dxm) * (1 - dym);
        ETAPSUM(i, j + 1) = ETAPSUM(i, j + 1) + etamm * wtmij1;
        PORSUM(i, j + 1) = PORSUM(i, j + 1) + porm(m) * wtmij1;
        SXXSUM(i, j + 1) = SXXSUM(i, j + 1) + sxxm(m) * wtmij1;
        SYYSUM(i, j + 1) = SYYSUM(i, j + 1) + syym(m) * wtmij1;
        GGGPSUM(i, j + 1) = GGGPSUM(i, j + 1) + 1 / gm(m) * wtmij1;
        ETAP0SUM(i, j + 1) = ETAP0SUM(i, j + 1) + etamm0 * wtmij1;
        ETAB0SUM(i, j + 1) = ETAB0SUM(i, j + 1) + etasmm0 * wtmij1;
        WTPSUM(i, j + 1) = WTPSUM(i, j + 1) + wtmij1;
        
        wtmi1j1 = (dxm) * (dym);
        ETAPSUM(i + 1, j + 1) = ETAPSUM(i + 1, j + 1) + etamm * wtmi1j1;
        PORSUM(i + 1, j + 1) = PORSUM(i + 1, j + 1) + porm(m) * wtmi1j1;
        SXXSUM(i + 1, j + 1) = SXXSUM(i + 1, j + 1) + sxxm(m) * wtmi1j1;
        SYYSUM(i + 1, j + 1) = SYYSUM(i + 1, j + 1) + syym(m) * wtmi1j1;
        GGGPSUM(i + 1, j + 1) = GGGPSUM(i + 1, j + 1) + 1 / gm(m) * wtmi1j1;
        ETAP0SUM(i + 1, j + 1) = ETAP0SUM(i + 1, j + 1) + etamm0 * wtmi1j1;
        ETAB0SUM(i + 1, j + 1) = ETAB0SUM(i + 1, j + 1) + etasmm0 * wtmi1j1;
        WTPSUM(i + 1, j + 1) = WTPSUM(i + 1, j + 1) + wtmi1j1;
        
        % Interpolate to Vx nodes
        % [i, j] -------- [i, j + 1]
        %   |                |
        %   |    o m         |
        %   |                |
        % [i + 1, j] ------- [i + 1, j + 1]
        % Indexes and distances
        j = fix((xm(m)) / dx) + 1;
        i = fix((ym(m) + dy / 2) / dy) + 1;
        if(j < 1)
            j = 1;
        elseif(j > Nx - 1)
            j = Nx - 1;
        end
        if(i < 1)
            i = 1;
        elseif(i > Ny)
            i = Ny;
        end
        dxm = (xm(m) - xvx(j)) / dx;
        dym = (ym(m) - yvx(i)) / dy;
        
        
        wtmij = (1 - dxm) * (1 - dym);
        RHOXSUM(i, j) = RHOXSUM(i, j) + rhomm * wtmij;
        ETADXSUM(i, j) = ETADXSUM(i, j) + 1 / (etadm) * wtmij;
        RHOFXSUM(i, j) = RHOFXSUM(i, j) + rhofm(m) * wtmij;
        VX0SUM(i, j) = VX0SUM(i, j) + vx0m(m) * rhomm * wtmij;
        PORXSUM(i, j) = PORXSUM(i, j) + porm(m) * wtmij;
        WTXSUM(i, j) = WTXSUM(i, j) + wtmij;
        
        
        wtmi1j = (1 - dxm) * (dym);
        RHOXSUM(i + 1, j) = RHOXSUM(i + 1, j) + rhomm * wtmi1j;
        ETADXSUM(i + 1, j) = ETADXSUM(i + 1, j) + 1 / (etadm) * wtmi1j;
        RHOFXSUM(i + 1, j) = RHOFXSUM(i + 1, j) + rhofm(m) * wtmi1j;
        VX0SUM(i + 1, j) = VX0SUM(i + 1, j) + vx0m(m) * rhomm * wtmi1j;
        PORXSUM(i + 1, j) = PORXSUM(i + 1, j) + porm(m) * wtmi1j;
        WTXSUM(i + 1, j) = WTXSUM(i + 1, j) + wtmi1j;
        
        
        wtmij1 = (dxm) * (1 - dym);
        RHOXSUM(i, j + 1) = RHOXSUM(i, j + 1) + rhomm * wtmij1;
        ETADXSUM(i, j + 1) = ETADXSUM(i, j + 1) + 1 / (etadm) * wtmij1;
        RHOFXSUM(i, j + 1) = RHOFXSUM(i, j + 1) + rhofm(m) * wtmij1;
        VX0SUM(i, j + 1) = VX0SUM(i, j + 1) + vx0m(m) * rhomm * wtmij1;
        PORXSUM(i, j + 1) = PORXSUM(i, j + 1) + porm(m) * wtmij1;
        WTXSUM(i, j + 1) = WTXSUM(i, j + 1) + wtmij1;
        
        
        wtmi1j1 = (dxm) * (dym);
        RHOXSUM(i + 1, j + 1) = RHOXSUM(i + 1, j + 1) + rhomm * wtmi1j1;
        ETADXSUM(i + 1, j + 1) = ETADXSUM(i + 1, j + 1) + 1 / (etadm) * wtmi1j1;
        RHOFXSUM(i + 1, j + 1) = RHOFXSUM(i + 1, j + 1) + rhofm(m) * wtmi1j1;
        VX0SUM(i + 1, j + 1) = VX0SUM(i + 1, j + 1) + vx0m(m) * rhomm * wtmi1j1;
        PORXSUM(i + 1, j + 1) = PORXSUM(i + 1, j + 1) + porm(m) * wtmi1j1;
        WTXSUM(i + 1, j + 1) = WTXSUM(i + 1, j + 1) + wtmi1j1;
        
        % Interpolate to Vy nodes
        % [i, j] -------- [i, j + 1]
        %   |                |
        %   |    o m         |
        %   |                |
        % [i + 1, j] ------- [i + 1, j + 1]
        % Indexes and distances
        j = fix((xm(m) + dx / 2) / dx) + 1;
        i = fix((ym(m)) / dy) + 1;
        if(j < 1)
            j = 1;
        elseif(j > Nx)
            j = Nx;
        end
        if(i < 1)
            i = 1;
        elseif(i > Ny - 1)
            i = Ny - 1;
        end
        dxm = (xm(m) - xvy(j)) / dx;
        dym = (ym(m) - yvy(i)) / dy;
        
        
        wtmij = (1 - dxm) * (1 - dym);
        RHOYSUM(i, j) = RHOYSUM(i, j) + rhomm * wtmij;
        ETADYSUM(i, j) = ETADYSUM(i, j) + 1 / (etadm) * wtmij;
        RHOFYSUM(i, j) = RHOFYSUM(i, j) + rhofm(m) * wtmij;
        VY0SUM(i, j) = VY0SUM(i, j) + vy0m(m) * rhomm * wtmij;
        PORYSUM(i, j) = PORYSUM(i, j) + porm(m) * wtmij;
        WTYSUM(i, j) = WTYSUM(i, j) + wtmij;
        
        
        wtmi1j = (1 - dxm) * (dym);
        RHOYSUM(i + 1, j) = RHOYSUM(i + 1, j) + rhomm * wtmi1j;
        ETADYSUM(i + 1, j) = ETADYSUM(i + 1, j) + 1 / (etadm) * wtmi1j;
        RHOFYSUM(i + 1, j) = RHOFYSUM(i + 1, j) + rhofm(m) * wtmi1j;
        VY0SUM(i + 1, j) = VY0SUM(i + 1, j) + vy0m(m) * rhomm * wtmi1j;
        PORYSUM(i + 1, j) = PORYSUM(i + 1, j) + porm(m) * wtmi1j;
        WTYSUM(i + 1, j) = WTYSUM(i + 1, j) + wtmi1j;
        
        
        wtmij1 = (dxm) * (1 - dym);
        RHOYSUM(i, j + 1) = RHOYSUM(i, j + 1) + rhomm * wtmij1;
        ETADYSUM(i, j + 1) = ETADYSUM(i, j + 1) + 1 / (etadm) * wtmij1;
        RHOFYSUM(i, j + 1) = RHOFYSUM(i, j + 1) + rhofm(m) * wtmij1;
        VY0SUM(i, j + 1) = VY0SUM(i, j + 1) + vy0m(m) * rhomm * wtmij1;
        PORYSUM(i, j + 1) = PORYSUM(i, j + 1) + porm(m) * wtmij1;
        WTYSUM(i, j + 1) = WTYSUM(i, j + 1) + wtmij1;
        
        
        wtmi1j1 = (dxm) * (dym);
        RHOYSUM(i + 1, j + 1) = RHOYSUM(i + 1, j + 1) + rhomm * wtmi1j1;
        ETADYSUM(i + 1, j + 1) = ETADYSUM(i + 1, j + 1) + 1 / (etadm) * wtmi1j1;
        RHOFYSUM(i + 1, j + 1) = RHOFYSUM(i + 1, j + 1) + rhofm(m) * wtmi1j1;
        VY0SUM(i + 1, j + 1) = VY0SUM(i + 1, j + 1) + vy0m(m) * rhomm * wtmi1j1;
        PORYSUM(i + 1, j + 1) = PORYSUM(i + 1, j + 1) + porm(m) * wtmi1j1;
        WTYSUM(i + 1, j + 1) = WTYSUM(i + 1, j + 1) + wtmi1j1;
        
    end % ends loop through markers
    
    % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    % end measuring marker time
    total_time_markers = total_time_markers + toc(t_marker);
    
    % start measuring nodes time
    t_nodes = tic;
    % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    % /////////////////////////////////////////////////////////////////////////////////////// 
    % third layer loop (timestep -> Ny -> Nx),  costly
    % /////////////////////////////////////////////////////////////////////////////////////// 
    
    % Computing ETA and RHO
    for i = 1:1:Ny
        % Basic nodes
        for j = 1:1:Nx
            if(WTSUM(i, j) > 0)
                RHO(i, j) = RHOSUM(i, j) / WTSUM(i, j);
                KKK(i, j) = (KKKSUM(i, j) / WTSUM(i, j));
                TTT(i, j) = (TTTSUM(i, j) / WTSUM(i, j));
                GGG(i, j) = 1 / (GGGSUM(i, j) / WTSUM(i, j));
                ETA0(i, j) = ETA0SUM(i, j) / WTSUM(i, j);
                COHT(i, j) = (COHTSUM(i, j) / WTSUM(i, j));
                FRIT(i, j) = (FRITSUM(i, j) / WTSUM(i, j));
                COHC(i, j) = (COHCSUM(i, j) / WTSUM(i, j));
                FRIC(i, j) = (FRICSUM(i, j) / WTSUM(i, j));
                DILC(i, j) = (DILCSUM(i, j) / WTSUM(i, j));
            end
        end
        % Vy nodes
        for j = 1:1:Nx1
            if(WTYSUM(i, j) > 0)
                RHOY(i, j) = RHOYSUM(i, j) / WTYSUM(i, j);
                ETADY(i, j) = 1 / (ETADYSUM(i, j) / WTYSUM(i, j));
                RHOFY(i, j) = RHOFYSUM(i, j) / WTYSUM(i, j);
                PORY(i, j) = PORYSUM(i, j) / WTYSUM(i, j);
            end
        end
    end
    
    for i = 1:1:Ny1
        % Vx nodes
        for j = 1:1:Nx
            if(WTXSUM(i, j) > 0)
                RHOX(i, j) = RHOXSUM(i, j) / WTXSUM(i, j);
                ETADX(i, j) = 1 / (ETADXSUM(i, j) / WTXSUM(i, j));
                RHOFX(i, j) = RHOFXSUM(i, j) / WTXSUM(i, j);
                PORX(i, j) = PORXSUM(i, j) / WTXSUM(i, j);
            end
        end
        %Pressure nodes
        for j = 1:1:Nx1
            if(WTPSUM(i, j) > 0)
                ETAP(i, j) = ETAPSUM(i, j) / WTPSUM(i, j);
                POR(i, j) = PORSUM(i, j) / WTPSUM(i, j);
                ETAB0(i, j) = ETAB0SUM(i, j) / WTPSUM(i, j) / POR(i, j);
                ETAP0(i, j) = ETAP0SUM(i, j) / WTPSUM(i, j);
                GGGP(i, j) = 1 / (GGGPSUM(i, j) / WTPSUM(i, j));
            end
        end
    end
    
    % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    % end measuring marker time
    total_nodes = total_nodes + toc(t_nodes);
    % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    % Save viscosity
    if(timestep == 1) % ??? only done once?
        ETA1 = ETA0;
        ETA = ETA0;
        ETA50 = ETA0;
    else
        ETA1 = ETA00; % variable ETA00 not defined until later
        ETA = ETA00;
        ETA50 = ETA00;
    end
    YNY00 = YNY0;
    
    OM = OM0; % ???
    
    % Multiple solving of equations
    if(yndtdecrease == 0)
        dt = max(min(dt * dtkoefup, dtelastic0), dtmin);
    else
        dt = max(min(dt, dtelastic0), dtmin);
    end
    
    yndtdecrease = 0;
    dt00 = dt;
    DSYLSQ = 0;
    ynlast = 0;
    
    % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    % start measuring global iteration time
    t_global_iterations = tic;
    % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    % /////////////////////////////////////////////////////////////////////////////////////// 
    % second layer loop (timestep -> number of global iterations),  maybe costly
    % /////////////////////////////////////////////////////////////////////////////////////// 
    
    for iterstep = 1:1:niterglobal
        
        % Limiting viscosity
        etamincur = shearmod * dt * 1e-4;
        
        % External P - nodes: symmetry
        pt(:, 1) = pt(:, 2);
        pt(:, Nx1) = pt(:, Nx);
        pt(1, :) = pt(2, :);
        pt(Ny1, :) = pt(Ny, :);
        pf(:, 1) = pf(:, 2);
        pf(:, Nx1) = pf(:, Nx);
        pf(1, :) = pf(2, :);
        pf(Ny1, :) = pf(Ny, :);
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % start measuring global iteration time
        t_plastic_strain = tic;
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % /////////////////////////////////////////////////////////////////////////////////////// 
        % fourth layer loop (timestep -> number of global iterations -> Ny -> Nx), 
        % very costly!!! reduce every single operation possible in here
        % /////////////////////////////////////////////////////////////////////////////////////// 
        
        % Basic nodes
        for i = 1:1:Ny
            for j = 1:1:Nx
                if(ETA(i, j) < etamincur)
                    ETA(i, j) = etamincur;
                end
                % Compute plastic strain rate
                if(ETA(i, j) < ETA0(i, j))
                    SIIB = (SXY(i, j)^2 + 0.5 * ((SXX(i, j) + SXX(i + 1, j) + SXX(i, j + 1) + SXX(i + 1, j + 1)) / 4)^2 + ...
                        0.5 * ((SYY(i, j) + SYY(i + 1, j) + SYY(i, j + 1) + SYY(i + 1, j + 1)) / 4)^2 + ...
                        0.5 * ((-SXX(i, j) - SYY(i, j) - SXX(i + 1, j) - SYY(i + 1, j)...
                         - SXX(i, j + 1) - SYY(i, j + 1) - SXX(i + 1, j + 1) - SYY(i + 1, j + 1)) / 4)^2)^0.5;
                    EIIB(i, j) = dy / faultwidth * (SIIB / 2 / ETA(i, j) - SIIB / 2 / ETA0(i, j));
                    IETAPLB(i, j) = (1 / ETA(i, j) - 1 / ETA0(i, j));
                else
                    EIIB(i, j) = 0;
                    IETAPLB(i, j) = 0;
                end
            end
        end
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % end measuring plastic strain time
        total_plastic_strain = total_plastic_strain + toc(t_plastic_strain);
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % start measuring global iteration time
        t_viscoplastic = tic;
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % /////////////////////////////////////////////////////////////////////////////////////// 
        % similar loop to above,  maybe merge them
        % fourth layer loop (timestep -> number of global iterations -> Ny -> Nx), 
        % very costly!!! reduce every single operation possible in here
        % /////////////////////////////////////////////////////////////////////////////////////// 
        
        % Computing viscosity and dilatation in pressure nodes
        for i = 2:1:Ny
            for j = 2:1:Nx
                % Compute viscoplastic viscosity
                IETAPL = (IETAPLB(i - 1, j - 1) + IETAPLB(i, j - 1) + IETAPLB(i - 1, j) + IETAPLB(i, j)) / 4;
                if(YNY0(i - 1, j - 1) > 0 || YNY0(i, j - 1) > 0 || YNY0(i - 1, j) > 0 || YNY0(i, j) > 0)
                    ETAP(i, j) = 1 / (1 / ETAP0(i, j) + IETAPL);
                    ETAB(i, j) = 1 / (1 / ETAB0(i, j) + dy / faultwidth * IETAPL * POR(i, j));
                else
                    ETAP(i, j) = ETAP0(i, j);
                    ETAB(i, j) = ETAB0(i, j);
                end
                % Check viscosity
                if(ETAP(i, j) < etamincur)
                    ETAP(i, j) = etamincur;
                end
                if(ETAB(i, j) * POR(i, j) < etamincur)
                    ETAB(i, j) = etamincur / POR(i, j);
                end
                % Pores compressibility
                GGGB(i, j) = GGGP(i, j) / POR(i, j);
                % Dilation
                % Zhao and Cai, International Journal of Rock Mechanics & Mining Sciences 47 (2010) 368–384
                % Weak sandstone parameters
                
                % /////////////////////////////////////////////////////////////////////////////////////// 
                % move variable declaration outside of the loop
                % /////////////////////////////////////////////////////////////////////////////////////// 
                
                ss3 = min(max((pt(i, j) - pf(i, j)) * 1e-6, 0), 100); % SIGMA3,  MPa
                aa = aa1 + aa2 * exp(-ss3 / aa3);
                bb = bb1 + bb2 * exp(-ss3 / bb3);
                cc = cc1 + cc2 / 100 * ss3^cc3;
                dilij = sin(aa * bb * (exp(-bb * gammapij) - exp(-cc * gammapij)) / (cc - bb) / 180 * pi);
                dili1j = sin(aa * bb * (exp(-bb * gammapi1j) - exp(-cc * gammapi1j)) / (cc - bb) / 180 * pi);
                dilij1 = sin(aa * bb * (exp(-bb * gammapij1) - exp(-cc * gammapij1)) / (cc - bb) / 180 * pi);
                dili1j1 = sin(aa * bb * (exp(-bb * gammapi1j1) - exp(-cc * gammapi1j1)) / (cc - bb) / 180 * pi);
                DILP(i, j) = 0; %2 * (dili1j1 * EIIB(i - 1, j - 1) + dilij1 * EIIB(i, j - 1) + dili1j * EIIB(i - 1, j) + dilij * EIIB(i, j)) / 4;
            end
        end
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % end measuring viscoplastic time
        total_viscoplastic = total_viscoplastic + toc(t_viscoplastic);
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        if(dt > 2e+5)
            pfscale = min(min(ETADX(2:Ny, 1:Nx))) * dx * 1e+17 / dt^2;
        else
            pfscale = min(min(ETADX(2:Ny, 1:Nx))) * dx;
        end
        
        ptscale = pfscale;
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % start measuring global matrices time
        t_global_matrices = tic;
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % /////////////////////////////////////////////////////////////////////////////////////// 
        % fourth layer loop (timestep -> number of global iterations -> Ny -> Nx), 
        % very costly!!! reduce every single operation possible in here
        % /////////////////////////////////////////////////////////////////////////////////////// 
        
        % 5)Composing global matrixes L(),  R()
        for j = 1:1:Nx1
            for i = 1:1:Ny1
                % Computing global indexes for vx, vy, p
                kp = ((j - 1) * Ny1 + (i - 1)) * 6 + 1;
                kx = kp + 1;
                ky = kp + 2;
                kpf = kp + 3;
                kxf = kp + 4;
                kyf = kp + 5;
                
                % 5a) Composing equation for vxs
                if(i == 1 || i == Ny1 || j == 1 || j == Nx || j == Nx1)
                    % Ghost nodes: 1 * vxs = 0
                    if(j == Nx1)
                        L(kx, kx) = 1;
                        R(kx) = 0;
                    end
                    
                    % Upper boundary
                    % prescribed velocity
                    if(i == 1 && j < Nx1)
                        L(kx, kx) = 1;
                        L(kx, kx + 6) = 1;
                        R(kx) = 2 * bcupper;
                    end
                    
                    % Lower boundary
                    % prescribed velocity
                    if(i == Ny1 && j < Nx1)
                        L(kx, kx) = 1;
                        L(kx, kx - 6) = 1;
                        R(kx) = 2 * bclower;
                    end
                    
                    % Left boundary:
                    if(j == 1 && i > 1 && i < Ny1)
                        
                        L(kx, kx) = 1;
                        L(kx, kx + 6 * Ny1) = -1;
                        R(kx) = 0;
                    end
                    
                    % Right boundary
                    if(j == Nx && i > 1 && i < Ny1)
                        L(kx, kx) = 1;
                        L(kx, kx - 6 * Ny1) = -1;
                        R(kx) = 0;
                    end
                    
                else
                    % Total X - Stokes: dSIGMAxxt' / dx + dSIGMAxyt' / dy - dPt / dx = -RHOt * gx
                    % SIGMAijt = 2 * ETA * EPSILONijs * K + SIGMAijt0 * (1 - K)
                    %             vxs2
                    %        vys1  |    vys3
                    %              |
                    %  vxs1 -- Pt1 -- vxs3 -- Pt2 -- vxs5
                    %              |
                    %        vys2  |    vys4
                    %             vxs4
                    % Viscosity
                    ETAXY1 = ETA(i - 1, j);
                    ETAXY2 = ETA(i, j);
                    ETAXX1 = ETAP(i, j);
                    ETAXX2 = ETAP(i, j + 1);
                    % Shear modulus
                    GXY1 = GGG(i - 1, j);
                    GXY2 = GGG(i, j);
                    GXX1 = GGGP(i, j);
                    GXX2 = GGGP(i, j + 1);
                    % Viscoelasticity factor
                    KXY1 = dt * GXY1 / (dt * GXY1 + ETAXY1);
                    KXY2 = dt * GXY2 / (dt * GXY2 + ETAXY2);
                    KXX1 = dt * GXX1 / (dt * GXX1 + ETAXX1);
                    KXX2 = dt * GXX2 / (dt * GXX2 + ETAXX2);
                    % Numerical viscosity
                    ETAXY1 = ETAXY1 * KXY1;
                    ETAXY2 = ETAXY2 * KXY2;
                    ETAXX1 = ETAXX1 * KXX1;
                    ETAXX2 = ETAXX2 * KXX2;
                    % Numerical stresses
                    SXY1 = SXY0(i - 1, j) * (1 - KXY1);
                    SXY2 = SXY0(i, j) * (1 - KXY2);
                    SXX1 = SXX0(i, j) * (1 - KXX1);
                    SXX2 = SXX0(i, j + 1) * (1 - KXX2);
                    % Density derivatives
                    dRHOdx = (RHOX(i, j + 1) - RHOX(i, j - 1)) / 2 / dx;
                    dRHOdy = (RHO(i, j) - RHO(i - 1, j)) / dy;
                    % Left part
                    L(kx, kx) = -4 / 3 * (ETAXX1 + ETAXX2) / dx^2 ...
                         - (ETAXY1 + ETAXY2) / dy^2 - gx * dt * dRHOdx - ascale * RHOX(i, j) / dt; %vxs3
                    L(kx, kx - Ny1 * 6) = 4 / 3 * ETAXX1 / dx^2; %vxs1
                    L(kx, kx + Ny1 * 6) = 4 / 3 * ETAXX2 / dx^2; %vxs5
                    L(kx, kx - 6) = ETAXY1 / dy^2; %vxs2
                    L(kx, kx + 6) = ETAXY2 / dy^2; %vxs4
                    L(kx, ky - 6) = ETAXY1 / dx / dy - 2 / 3 * ETAXX1 / dx / dy - gx * dt * dRHOdy / 4; %vys1
                    L(kx, ky) = -ETAXY2 / dx / dy + 2 / 3 * ETAXX1 / dx / dy - gx * dt * dRHOdy / 4; %vys2
                    L(kx, ky - 6 + Ny1 * 6) = -ETAXY1 / dx / dy + 2 / 3 * ETAXX2 / dx / dy - gx * dt * dRHOdy / 4; %vys3
                    L(kx, ky + Ny1 * 6) = ETAXY2 / dx / dy - 2 / 3 * ETAXX2 / dx / dy - gx * dt * dRHOdy / 4; %vys4
                    L(kx, kp) = ptscale / dx; %Pt1'
                    L(kx, kp + Ny1 * 6) = -ptscale / dx; %Pt2'
                    % Right part
                    R(kx) = -RHOX(i, j) * (ascale * VX0(i, j) / dt + gx) - (SXX2 - SXX1) / dx - (SXY2 - SXY1) / dy;
                end
                
                % 5b) Composing equation for vys
                if(j == 1 || j == Nx1 || i == 1 || i == Ny || i == Ny1)
                    % Ghost nodes: 1 * vys = 0
                    if(i == Ny1)
                        L(ky, ky) = 1;
                        R(ky) = 0;
                    end
                    
                    % Left boundary
                    % Free Slip
                    if(j == 1)
                        L(ky, ky) = 1;
                        L(ky, ky + Ny1 * 6) = 1;
                        R(ky) = 0;
                    end
                    
                    % Right boundary
                    % Free Slip
                    if(j == Nx1)
                        L(ky, ky) = 1;
                        L(ky, ky - Ny1 * 6) = 1;
                        R(ky) = 0;
                    end
                    
                    % Upper boundary: no penetration
                    if(i == 1 && j > 1 && j < Nx1)
                        L(ky, ky) = 1;
                        R(ky) = 0;
                    end
                    
                    % Lower boundary: no penetration
                    if(i == Ny && j > 1 && j < Nx1)
                        L(ky, ky) = 1;
                        R(ky) = 0;
                    end
                    
                else
                    % Total Y - Stokes: dSIGMAyxt' / dx + dSIGMAyyt' / dy - dPt / dy = -RHOt * gy
                    % y - Stokes equation: dSIGMA'yx / dx + dSIGMA'yy / dy - dP / dy = -RHO * gy
                    %
                    %               vys2
                    %                |
                    %         vxs1  Pt1  vxs3
                    %                |
                    %   vys1 --------- vys3 -------- vys5
                    %                |
                    %         vxs2  Pt2  vxs4
                    %                |
                    %               vys4
                    % Viscosity
                    ETAXY1 = ETA(i, j - 1);
                    ETAXY2 = ETA(i, j);
                    ETAYY1 = ETAP(i, j);
                    ETAYY2 = ETAP(i + 1, j);
                    % Shear modulus
                    GXY1 = GGG(i, j - 1);
                    GXY2 = GGG(i, j);
                    GYY1 = GGGP(i, j);
                    GYY2 = GGGP(i + 1, j);
                    % Viscoelasticity factor
                    KXY1 = dt * GXY1 / (dt * GXY1 + ETAXY1);
                    KXY2 = dt * GXY2 / (dt * GXY2 + ETAXY2);
                    KYY1 = dt * GYY1 / (dt * GYY1 + ETAYY1);
                    KYY2 = dt * GYY2 / (dt * GYY2 + ETAYY2);
                    % Numerical viscosity
                    ETAXY1 = ETAXY1 * KXY1;
                    ETAXY2 = ETAXY2 * KXY2;
                    ETAYY1 = ETAYY1 * KYY1;
                    ETAYY2 = ETAYY2 * KYY2;
                    % Numerical stresses
                    SXY1 = SXY0(i, j - 1) * (1 - KXY1);
                    SXY2 = SXY0(i, j) * (1 - KXY2);
                    SYY1 = SYY0(i, j) * (1 - KYY1);
                    SYY2 = SYY0(i + 1, j) * (1 - KYY2);
                    % Density derivatives
                    dRHOdy = (RHOY(i + 1, j) - RHOY(i - 1, j)) / 2 / dy;
                    dRHOdx = (RHO(i, j) - RHO(i, j - 1)) / dx;
                    % Left part
                    L(ky, ky) = -4 / 3 * (ETAYY1 + ETAYY2) / dy^2 - ...
                        (ETAXY1 + ETAXY2) / dx^2 - gy * dt * dRHOdy - ascale * RHOY(i, j) / dt; %vys3
                    L(ky, ky - Ny1 * 6) = ETAXY1 / dx^2; %vys1
                    L(ky, ky + Ny1 * 6) = ETAXY2 / dx^2; %vys5
                    L(ky, ky - 6) = 4 / 3 * ETAYY1 / dy^2; %vys2
                    L(ky, ky + 6) = 4 / 3 * ETAYY2 / dy^2; %vys4
                    L(ky, kx - Ny1 * 6) = ETAXY1 / dx / dy - 2 / 3 * ETAYY1 / dx / dy - gy * dt * dRHOdx / 4; %vxs1
                    L(ky, kx + 6 - Ny1 * 6) = -ETAXY1 / dx / dy + 2 / 3 * ETAYY2 / dx / dy - gy * dt * dRHOdx / 4; %vxs2
                    L(ky, kx) = -ETAXY2 / dx / dy + 2 / 3 * ETAYY1 / dx / dy - gy * dt * dRHOdx / 4; %vxs3
                    L(ky, kx + 6) = ETAXY2 / dx / dy - 2 / 3 * ETAYY2 / dx / dy - gy * dt * dRHOdx / 4; %vxs4
                    L(ky, kp) = ptscale / dy; %Pt1'
                    L(ky, kp + 6) = -ptscale / dy; %Pt2'
                    % Right part
                    R(ky) = -RHOY(i, j) * (ascale * VY0(i, j) / dt + gy) - (SYY2 - SYY1) / dy - (SXY2 - SXY1) / dx;
                end
                
                
                % 5c) Composing equation for Pt
                if(i == 1 || j == 1 || i == Ny1 || j == Nx1)% || (i == 2 && j == 2))
                    % BC equation: 1 * Pt = 0
                    L(kp, kp) = 1;
                    R(kp) = 0;
                else
                    % Solid Continuity: dVxs / dx + dVys / dy + (Pt - Pf) / ETAbulk = 0
                    %              vys1
                    %               |
                    %        vxs1 -- Pt, Pf -- vxs2
                    %               |
                    %              vys2
                    % Drained compressibility
                    BETTADRAINED = (1 / GGGB(i, j) + BETTASOLID) / (1 - POR(i, j));
                    % Biott - Willis koefficient
                    KBW = 1 - BETTASOLID / BETTADRAINED;
                    % Left part
                    L(kp, kx - Ny1 * 6) = -1 / dx; %vxs1
                    L(kp, kx) = 1 / dx; %vxs2
                    L(kp, ky - 6) = -1 / dy; %vys1
                    L(kp, ky) = 1 / dy; %vys2
                    L(kp, kp) = ptscale * (1 / ETAB(i, j) / (1 - POR(i, j)) + gggbkoef * BETTADRAINED / dt); %Pt
                    L(kp, kpf) = -pfscale * (1 / ETAB(i, j) / (1 - POR(i, j)) + gggbkoef * BETTADRAINED * KBW / dt); %Pf
                    % Right part
                    R(kp) = gggbkoef * BETTADRAINED * (PT0(i, j) - KBW * PF0(i, j)) / dt + DILP(i, j);
                end
                
                % 5d) Composing equation for vxD
                if(i == 1 || i == Ny1 || j == 1 || j == Nx || j == Nx1)
                    % Ghost nodes: 1 * vxs = 0
                    if(j == Nx1)
                        L(kxf, kxf) = 1;
                        R(kxf) = 0;
                    end
                    
                    % Upper boundary: symmetry
                    if(i == 1 && j < Nx1)
                        L(kxf, kxf) = 1;
                        L(kxf, kxf + 6) = -1;
                        R(kxf) = 0;
                    end
                    
                    % Lower boundary: symmetry
                    if(i == Ny1 && j < Nx1)
                        L(kxf, kxf) = 1;
                        L(kxf, kxf - 6) = -1;
                        R(kxf) = 0;
                    end
                    
                    % Left boundary
                    % no penetration
                    if(j == 1)
                        L(kxf, kxf) = 1;
                        R(kxf) = 0;%bcvxfleft;
                    end
                    
                    % Right boundary
                    % no penetration
                    if(j == Nx)
                        L(kxf, kxf) = 1;
                        R(kxf) = 0;
                    end
                    
                else
                    % Fluid X - Darsi:  - ETAfluid / K * VxD - dPf / dx = -RHOf * gx + RHOf * DVxs / Dt
                    %
                    %  Pf1 --- vxD, vxs --- Pf2
                    %
                    % Left part
                    L(kxf, kxf) = -ETADX(i, j) - RHOFX(i, j) / PORX(i, j) * ascale / dt; %vxD
                    L(kxf, kx) = -RHOFX(i, j) * ascale / dt; %vxs
                    L(kxf, kpf) = pfscale / dx; %Pf1'
                    L(kxf, kpf + Ny1 * 6) = -pfscale / dx; %Pf2'
                    % Right part
                    R(kxf) = -RHOFX(i, j) * (ascale * VXF0(i, j) / dt + gx);
                end
                
                % 5e) Composing equation for vyD
                if(j == 1 || j == Nx1 || i == 1 || i == Ny || i == Ny1)
                    % Ghost nodes: 1 * vxs = 0
                    if(i == Ny1)
                        L(kyf, kyf) = 1;
                        R(kyf) = 0;
                    end
                    
                    % Left boundary
                    % symmetry
                    if(j == 1 && i > 1 && i < Ny)
                        L(kyf, kyf) = 1;
                        L(kyf, kyf + Ny1 * 6) = -1;
                        R(kyf) = 0;
                    end
                    
                    % Right boundary
                    % symmetry
                    if(j == Nx1 && i > 1 && i < Ny)
                        L(kyf, kyf) = 1;
                        L(kyf, kyf - Ny1 * 6) = -1;
                        R(kyf) = 0;
                    end
                    
                    % Upper boundary: no penetration
                    if(i == 1)
                        L(kyf, kyf) = 1;
                        R(kyf) = bcvyflower;
                    end
                    
                    % Lower boundary: no penetration
                    if(i == Ny)
                        L(kyf, kyf) = 1;
                        R(kyf) = bcvyflower;
                    end
                else
                    % Fluid Y - Darsi:  - ETAfluid / K * VyD - dPf / dy = -RHOf * gy + RHOf * DVys / Dt
                    %
                    %   Pf1
                    %    |
                    %   vyD, vy
                    %    |
                    %   Pf2
                    %
                    % Left part
                    L(kyf, kyf) = -ETADY(i, j) - RHOFY(i, j) / PORY(i, j) * ascale / dt; %vyD
                    L(kyf, ky) = -RHOFY(i, j) * ascale / dt; %vys
                    L(kyf, kpf) = pfscale / dy; %Pf1'
                    L(kyf, kpf + 6) = -pfscale / dy; %Pf2'
                    % Right part
                    R(kyf) = -RHOFY(i, j) * (ascale * VYF0(i, j) / dt + gy);
                end
                
                
                % 5f) Composing equation for Pf
                if(i == 1 || j == 1 || i == Ny1 || j == Nx1 || i == 2 || i == Ny)
                    % BC equation: 1 * Pf = 0
                    L(kpf, kpf) = 1;
                    R(kpf) = 0;
                    
                    % Real BC
                    if(i == 2 && (j > 2 || gggbkoef == 1)) % gggbkoef always defined as 1
                        L(kpf, kpf) = 1 * pfscale;
                        L(kpf, kp) = -1 * ptscale;
                        R(kpf) = -PTFDIFF;
                    end
                    if(i == 2 && j == 2 && gggbkoef == 0) % gggbkoef always defined as 1
                        L(kpf, kpf) = 1 * pfscale;
                        L(kpf, kp) = 0;
                        R(kpf) = PCONF;
                    end
                    
                    % Real BC
                    if(i == Ny && (j < Nx || gggbkoef == 1)) % gggbkoef always defined as 1
                        L(kpf, kpf) = 1 * pfscale;
                        L(kpf, kp) = -1 * ptscale;
                        R(kpf) = -PTFDIFF;
                    end
                    if(i == Ny && j == Nx && gggbkoef == 0) % gggbkoef always defined as 1
                        L(kpf, kpf) = 1 * pfscale;
                        L(kpf, kp) = 0;
                        R(kpf) = PCONF;
                    end
                    
                else
                    % Fluid Continuity: dVxD / dx + dVyD / dy - (Pt - Pf) / ETAbulk = 0
                    %              vyD1
                    %               |
                    %        vxD1 -- Pt, Pf -- vxD2
                    %               |
                    %              vyD2
                    % Compute elastic coefficients
                    % Drained compressibility
                    BETTADRAINED = (1 / GGGB(i, j) + BETTASOLID) / (1 - POR(i, j));
                    % Biott - Willis koefficient
                    KBW = 1 - BETTASOLID / BETTADRAINED;
                    % Skempton koefficient
                    KSK = (BETTADRAINED - BETTASOLID) / (BETTADRAINED - BETTASOLID + POR(i, j) * (BETTAFLUID - BETTASOLID));
                    % Left part
                    L(kpf, kxf - Ny1 * 6) = -1 / dx; %vxs1
                    L(kpf, kxf) = 1 / dx; %vxs2
                    L(kpf, kyf - 6) = -1 / dy; %vys1
                    L(kpf, kyf) = 1 / dy; %vys2
                    L(kpf, kp) = -ptscale * (1 / ETAB(i, j) / (1 - POR(i, j)) + BETTADRAINED * KBW / dt); %Pt
                    L(kpf, kpf) = pfscale * (1 / ETAB(i, j) / (1 - POR(i, j)) + BETTADRAINED * KBW / KSK / dt); %Pf
                    % Right part
                    R(kpf) = -BETTADRAINED * KBW * (PT0(i, j) - 1 / KSK * PF0(i, j)) / dt - DILP(i, j);
                end
            end
        end
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % end measuring global matrices time
        total_global_matrices = total_global_matrices + toc(t_global_matrices);
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % start measuring solve matrix time
        t_solve = tic;
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % 6) Solving matrix
        S = L\R; % This line causes about 60% of all computational cost
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % end measuring solve matrix time
        total_solve = total_solve+toc(t_solve);
        fprintf('%d \n', total_solve);
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % start measuring reload solution time
        t_reload_solution = tic;
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % /////////////////////////////////////////////////////////////////////////////////////// 
        % fourth layer loop (timestep -> number of global iterations -> Nx -> Ny), 
        % very costly!!! reduce every single operation possible in here
        % /////////////////////////////////////////////////////////////////////////////////////// 
        
        % 7) Reload solution
        % pfavr = 0;
        % pcount = 0;
        for j = 1:1:Nx1
            for i = 1:1:Ny1
                % Global indexes for vx, vy, P
                kp = ((j - 1) * Ny1 + (i - 1)) * 6 + 1; % start loops at 0  - >  kp = (j * Ny1 + i) * 6 + 1  = >  -2 operations each time
                kx = kp + 1;
                ky = kp + 2;
                kpf = kp + 3;
                kxf = kp + 4;
                kyf = kp + 5;
                % Reload solution
                pt(i, j) = S(kp) * ptscale;
                vxs(i, j) = S(kx);
                vys(i, j) = S(ky);
                pf(i, j) = S(kpf) * pfscale;
                vxD(i, j) = S(kxf);
                vyD(i, j) = S(kyf);
            end
        end
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % end measuring reload solution time
        total_reload_solution = total_reload_solution + toc(t_reload_solution);
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        Vmax = max(max(VSLIPB));
        
        if (dt > 1e4 && Vmax < 1e-7)
            pres_correction = 1;
            [length_pt,  width_pt] = size(pt);
            avgpt = sum(sum(pt)) / (length_pt * width_pt); %calculate average total pressure
            diffpt = (PCONF + PTFDIFF) - avgpt;
            pt(:, :) = pt(:, :) + diffpt;
        else
            pres_correction = 0;
        end
        
        % Velocity change
        DVX0 = vxs - VX0;
        DVY0 = vys - VY0;
        
        % Define timestep
        dt0 = dt;
        yn = 0;
        
        % /////////////////////////////////////////////////////////////////////////////////////// 
        % fdefine variables outside of loops and set zero before using
        % /////////////////////////////////////////////////////////////////////////////////////// 
        
        % Plastic iterations
        % Compute strain rate,  stress and stress change
        ESP = zeros(Ny, Nx);
        EXY = zeros(Ny, Nx);
        SXY = zeros(Ny, Nx);
        DSXY = zeros(Ny, Nx);
        EXX = zeros(Ny1, Nx1);
        SXX = zeros(Ny1, Nx1);
        DSXX = zeros(Ny1, Nx1);
        EYY = zeros(Ny1, Nx1);
        SYY = zeros(Ny1, Nx1);
        DSYY = zeros(Ny1, Nx1);
        EII = zeros(Ny1, Nx1);
        EIIVP = zeros(Ny1, Nx1);
        SII = zeros(Ny1, Nx1);
        DIS = zeros(Ny1, Nx1);
        
        EL_DECOM = zeros(Ny1, Nx1);   % Elastic (de)compaction
        VIS_COMP = zeros(Ny1, Nx1);   % Viscous compaction
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % start measuring internal nodes time
        t_internal_nodes = tic;
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % /////////////////////////////////////////////////////////////////////////////////////// 
        % fourth layer loop (timestep -> number of global iterations -> Ny -> Nx), 
        % very costly!!! reduce every single operation possible in here
        % /////////////////////////////////////////////////////////////////////////////////////// 
        
        % Process internal basic nodes
        for i = 1:1:Ny
            for j = 1:1:Nx
                % ESP = 1 / 2(dVy / dx - dVx / dy),  EXY,  SXY,  DSXY
                ESP(i, j) = 1 / 2 * ((vys(i, j + 1) - vys(i, j)) / dx - (vxs(i + 1, j) - vxs(i, j)) / dy); % move " / 2" outside of loop
                EXY(i, j) = 1 / 2 * ((vxs(i + 1, j) - vxs(i, j)) / dy + (vys(i, j + 1) - vys(i, j)) / dx); % move " / 2" outside of loop
                KXY = dt * GGG(i, j) / (dt * GGG(i, j) + ETA(i, j));
                SXY(i, j) = 2 * ETA(i, j) * EXY(i, j) * KXY + SXY0(i, j) * (1 - KXY); % move " * 2" outside of loop
                DSXY(i, j) = SXY(i, j) - SXY0(i, j);
            end
        end
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % end measuring internal nodes time
        total_internal_nodes = total_internal_nodes + toc(t_internal_nodes);
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % start measuring pressure cell time
        t_pressure_cell = tic;
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % /////////////////////////////////////////////////////////////////////////////////////// 
        % very similar loop as above,  maybe merge
        % fourth layer loop (timestep -> number of global iterations -> Ny -> Nx), 
        % very costly!!! reduce every single operation possible in here
        % /////////////////////////////////////////////////////////////////////////////////////// 
        
        % Process pressure cells
        for i = 2:1:Ny
            for j = 2:1:Nx
                % EXX,  SXX,  DSXX
                EXX(i, j) = (2 * (vxs(i, j) - vxs(i, j - 1)) / dx - (vys(i, j) - vys(i - 1, j)) / dy) / 3; % move " / 3" outside of loop
                EYY(i, j) = (2 * (vys(i, j) - vys(i - 1, j)) / dy - (vxs(i, j) - vxs(i, j - 1)) / dx) / 3; % move " / 3" outside of loop
                KXX = dt * GGGP(i, j) / (dt * GGGP(i, j) + ETAP(i, j));
                SXX(i, j) = 2 * ETAP(i, j) * EXX(i, j) * KXX + SXX0(i, j) * (1 - KXX);
                SYY(i, j) = 2 * ETAP(i, j) * EYY(i, j) * KXX + SYY0(i, j) * (1 - KXX);
                DSXX(i, j) = SXX(i, j) - SXX0(i, j);
                DSYY(i, j) = SYY(i, j) - SYY0(i, j);
            end
        end
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % end measuring pressure cell time
        total_pressure_cell = total_pressure_cell + toc(t_pressure_cell);
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % External P - nodes: symmetry
        pt(:, 1) = pt(:, 2);
        pt(:, Nx1) = pt(:, Nx);
        pt(1, :) = pt(2, :);
        pt(Ny1, :) = pt(Ny, :);
        pf(:, 1) = pf(:, 2);
        pf(:, Nx1) = pf(:, Nx);
        pf(1, :) = pf(2, :);
        pf(Ny1, :) = pf(Ny, :);
        EXX(:, 1) = EXX(:, 2);
        EXX(:, Nx1) = EXX(:, Nx);
        EXX(1, :) = EXX(2, :);
        EXX(Ny1, :) = EXX(Ny, :);
        SXX(:, 1) = SXX(:, 2);
        SXX(:, Nx1) = SXX(:, Nx);
        SXX(1, :) = SXX(2, :);
        SXX(Ny1, :) = SXX(Ny, :);
        SXX0(:, 1) = SXX0(:, 2);
        SXX0(:, Nx1) = SXX0(:, Nx);
        SXX0(1, :) = SXX0(2, :);
        SXX0(Ny1, :) = SXX0(Ny, :);
        EYY(:, 1) = EYY(:, 2);
        EYY(:, Nx1) = EYY(:, Nx);
        EYY(1, :) = EYY(2, :);
        EYY(Ny1, :) = EYY(Ny, :);
        SYY(:, 1) = SYY(:, 2);
        SYY(:, Nx1) = SYY(:, Nx);
        SYY(1, :) = SYY(2, :);
        SYY(Ny1, :) = SYY(Ny, :);
        SYY0(:, 1) = SYY0(:, 2);
        SYY0(:, Nx1) = SYY0(:, Nx);
        SYY0(1, :) = SYY0(2, :);
        SYY0(Ny1, :) = SYY0(Ny, :);
        ETAP(:, 1) = ETAP(:, 2);
        ETAP(:, Nx1) = ETAP(:, Nx);
        ETAP(1, :) = ETAP(2, :);
        ETAP(Ny1, :) = ETAP(Ny, :);
        ETAB(:, 1) = ETAB(:, 2);
        ETAB(:, Nx1) = ETAB(:, Nx);
        ETAB(1, :) = ETAB(2, :);
        ETAB(Ny1, :) = ETAB(Ny, :);
        GGGP(:, 1) = GGGP(:, 2);
        GGGP(:, Nx1) = GGGP(:, Nx);
        GGGP(1, :) = GGGP(2, :);
        GGGP(Ny1, :) = GGGP(Ny, :);
        GGGB(:, 1) = GGGB(:, 2);
        GGGB(:, Nx1) = GGGB(:, Nx);
        GGGB(1, :) = GGGB(2, :);
        GGGB(Ny1, :) = GGGB(Ny, :);
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % start measuring pressure cell time
        t_pressure_cell2 = tic;
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % /////////////////////////////////////////////////////////////////////////////////////// 
        % fourth layer loop (timestep -> number of global iterations -> Ny -> Nx), 
        % very costly!!! reduce every single operation possible in here
        % /////////////////////////////////////////////////////////////////////////////////////// 
        
        % Compute stress and strain rate invariants and dissipation
        % Process pressure cells
        for i = 2:1:Ny
            for j = 2:1:Nx
                % EXY term is averaged from four surrounding basic nodes
                EXY2 = (EXY(i, j)^2 + EXY(i - 1, j)^2 + EXY(i, j - 1)^2 + EXY(i - 1, j - 1)^2) / 4;
                EII(i, j) = (0.5 * (EXX(i, j)^2 + EYY(i, j)^2) + EXY2)^0.5;
                EXYVP2 = ((SXY(i, j) / 2 / ETA(i, j))^2 + (SXY(i - 1, j) / 2 / ETA(i - 1, j))^2 + (SXY(i, j - 1) / 2 / ETA(i, j - 1))^2 + (SXY(i - 1, j - 1) / 2 / ETA(i - 1, j - 1))^2) / 4;
                EIIVP(i, j) = (0.5 * ((SXX(i, j) / 2 / ETAP(i, j))^2 + (SYY(i, j) / 2 / ETAP(i, j))^2) + EXYVP2)^0.5;
                % Second strain rate invariant SII
                % SXY term is averaged from four surrounding basic nodes
                SXY2 = (SXY(i, j)^2 + SXY(i - 1, j)^2 + SXY(i, j - 1)^2 + SXY(i - 1, j - 1)^2) / 4;
                SII(i, j) = (0.5 * (SXX(i, j)^2 + SYY(i, j)^2) + SXY2)^0.5;
                
                % Dissipation
                DISXY = (SXY(i, j)^2 / 2 / ETA(i, j) + SXY(i - 1, j)^2 / 2 / ETA(i - 1, j) + SXY(i, j - 1)^2 / 2 / ETA(i, j - 1) + SXY(i - 1, j - 1)^2 / 2 / ETA(i - 1, j - 1)) / 4;
                DIS(i, j) = SXX(i, j)^2 / 2 / ETAP(i, j) + SYY(i, j)^2 / 2 / ETAP(i, j) + 2 * DISXY;
            end
        end
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % end measuring pressure cell time
        total_pressure_cell2 = total_pressure_cell2 + toc(t_pressure_cell2);
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % Update viscosity for yielding
        AXY = zeros(Ny, Nx);
        % dt0 = dt;dt = dt * 1.1;
        ETA5 = ETA0;
        % Basic nodes
        DSY = zeros(Ny, Nx);
        YNY = zeros(Ny, Nx);
        SigmaY = zeros(Ny, Nx);
        SII_fault = zeros(Ny, Nx);
        ynpl = 0;
        ddd = 0;
        dtrobert = dt;
        dtlapusta = 1e7;
        OM5 = OM;
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % start measuring global iteration time
        t_plasticity = tic;
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % /////////////////////////////////////////////////////////////////////////////////////// 
        % fourth layer loop (timestep -> number of global iterations -> Ny -> Nx), 
        % possibly very costly!!! reduce every single operation possible in here
        % /////////////////////////////////////////////////////////////////////////////////////// 
        
        % Power law plasticity model Yi et al.,  2018
        % Journal of Offshore Mechanics and Arctic Engineering
        dtslip = 1e+30;
        if(timestep > tyield)
            for i = 1:1:Ny
                if(y(i)  >= upper_block && y(i)  <= lower_block)
                    for j = 1:1:Nx
                        % SXX,  pt are averaged from four surrounding pressure nodes
                        SIIB(i, j) = (SXY(i, j)^2 + 1 / 2 * ((SXX(i, j) + SXX(i + 1, j) + SXX(i, j + 1) + SXX(i + 1, j + 1)) / 4)^2 + ...
                            1 / 2 * ((SYY(i, j) + SYY(i + 1, j) + SYY(i, j + 1) + SYY(i + 1, j + 1)) / 4)^2 + ...
                            1 / 2 * (( - SXX(i, j) - SYY(i, j) - SXX(i + 1, j) - SYY(i + 1, j)...
                             - SXX(i, j + 1) - SYY(i, j + 1) - SXX(i + 1, j + 1) - SYY(i + 1, j + 1)) / 4)^2)^0.5;
                        ptB = (pt(i, j) + pt(i + 1, j) + pt(i, j + 1) + pt(i + 1, j + 1)) / 4;
                        pfB = (pf(i, j) + pf(i + 1, j) + pf(i, j + 1) + pf(i + 1, j + 1)) / 4;
                        % Computing "elastic" stress invariant
                        kfxy = ETA(i, j) / (GGG(i, j) * dt + ETA(i, j));
                        siiel = SIIB(i, j) / kfxy;
                        % Compute maximal stress invariant with viscous viscosity
                        kfxy0 = ETA0(i, j) / (GGG(i, j) * dt + ETA0(i, j));
                        SIIB0 = siiel * kfxy0;
                        
                        % Assign PEFFB
                        PEFFB0(i, j) = (ptB - pfB);
                        PEFFB(i, j) = (ptB - pfB);
                        PEFFB1(i, j) = PEFFB(i, j) - PEFFB0(i, j);
                        % Compute old viscoplastic slip rate
                        % Compute PEFF
                        prB = (ptB - pfB);
                        
                        if(prB < 1e3)
                            prB = 1e3;
                        end
                        % Compute old power law strain rate
                        SIIB1 = SIIB(i, j);
                        
                        %Compute slip velocity for current stress invariant and state
                        V = 2 * V0 * sinh(max((SIIB1), 0) / ARSF(i, j) / prB) * ...
                            exp( - (BRSF(i, j) * OM(i, j) + FRIC(i, j)) / ARSF(i, j));
                        
                        EIISLIP = V / dx / 2;
                        
                        % Compute new ETAVP
                        ETAPL = SIIB1 / 2 / EIISLIP;
                        ETAVP = 1 / (1 / ETA0(i, j) + 1 / ETAPL);
                        % Compute new stress invariant
                        kfxy1 = ETAVP / (GGG(i, j) * dt + ETAVP);
                        SIIB2 = siiel * kfxy1;
                        DSIIB1 = SIIB2 - SIIB1;
                        
                        
                        %Compute slip velocity for current stress invariant and state
                        V = 2 * V0 * sinh(max((SIIB2), 0) / ARSF(i, j) / prB) * ...
                            exp( - (BRSF(i, j) * OM(i, j) + FRIC(i, j)) / ARSF(i, j));
                        
                        EIISLIP = V / dx / 2;
                        
                        % Compute new ETAVP
                        ETAPL = SIIB2 / 2 / EIISLIP;
                        ETAVP = 1 / (1 / ETA0(i, j) + 1 / ETAPL);
                        % Compute new stress invariant
                        kfxy1 = ETAVP / (GGG(i, j) * dt + ETAVP);
                        SIIB3 = siiel * kfxy1;
                        DSIIB2 = SIIB3 - SIIB2;
                        
                        if((DSIIB1 >= 0 && DSIIB2 <= 0) || (DSIIB1 <= 0 && DSIIB2 >= 0))
                            DSIIB = 1e+9;
                            ijk = 0;
                            
                            % /////////////////////////////////////////////////////////////////////////////////////// 
                            % fifth layer loop (timestep -> number of global iterations -> Ny -> Nx -> while(), 
                            % possibly extremely costly!!! reduce every single operation possible in here
                            % /////////////////////////////////////////////////////////////////////////////////////// 
                            
                            while(abs(DSIIB) > 1e-3)
                                SIIB4 = (SIIB1 + SIIB2) / 2;
                                
                                %Compute slip velocity for current stress invariant and state
                                V = 2 * V0 * sinh(max((SIIB4), 0) / ARSF(i, j) / prB) * ...
                                    exp( - (BRSF(i, j) * OM(i, j) + FRIC(i, j)) / ARSF(i, j));
                                
                                EIISLIP = V / dx / 2;
                                
                                % Compute new ETAVP
                                ETAPL = SIIB4 / 2 / EIISLIP;
                                ETAVP = 1 / (1 / ETA0(i, j) + 1 / ETAPL);
                                % Compute new stress invariant
                                kfxy1 = ETAVP / (GGG(i, j) * dt + ETAVP);
                                SIIB5 = siiel * kfxy1;
                                DSIIB = SIIB5 - SIIB4;
                                if((DSIIB >= 0 && DSIIB1 >= 0) || (DSIIB <= 0 && DSIIB1 <= 0))
                                    SIIB1 = SIIB4;
                                else
                                    SIIB2 = SIIB4;
                                end
                                ijk = ijk + 1;
                            end
                        end
                        
                        
                        if(V * dt / LRSF(i, j) > 1e-6)
                            OM5(i, j) = log(V0 / V + (exp(OM0(i, j)) - V0 / V) * exp( - V * dt / LRSF(i, j)));
                        else
                            OM5(i, j) = log(exp(OM0(i, j)) * (1 - V * dt / LRSF(i, j)) + V0 * dt / LRSF(i, j));
                        end
                        
                        
                        kfxy1 = ETAVP / (GGG(i, j) * dt + ETAVP);
                        SIGMA2 = siiel * kfxy1;
                        
                        % Compute yielding stress
                        syield = max(syieldmin, (ptB - pfB) * ARSF(i, j) * asinh(V / 2 / V0 * exp((BRSF(i, j) * OM5(i, j) + FRIC(i, j)) / ARSF(i, j))));
                        
                        % Compute visco - plastic viscosity
                        D = 1;
                        etapl = ETA0(i, j) * syield / (ETA0(i, j) * V / D + syield);
                        
                        %LDZ: save syield
                        SigmaY(i, j) = syield;
                        VSLIPB(i, j) = V;
                        SII_fault(i, j) = SIIB4;
                        
                        
                        % Timestep criterion,  Lapusta et al.,  2000; Lapusta and Liu,  2009
                        B = 1 / BETTASOLID;
                        vi = (3 * B - 2 * GGG(i, j)) / (6 * B + 2 * GGG(i, j));
                        k = 2 / pi * GGG(i, j) / (1 - vi) / dx;
                        xi = 1 / 4 * (k * LRSF(i, j) / ARSF(i, j) / prB - (BRSF(i, j) - ...
                            ARSF(i, j)) / ARSF(i, j))^2 - k * LRSF(i, j) / ARSF(i, j) / prB;
                        if(xi < 0)
                            dTETAmax = min(1 - (BRSF(i, j) - ARSF(i, j)) * prB / (k * LRSF(i, j)), 0.2);
                        else
                            dTETAmax = min(ARSF(i, j) * prB / (k * LRSF(i, j) - (BRSF(i, j) - ARSF(i, j)) * prB), 0.2);
                        end
                        dtlapusta = min(dtlapusta, dTETAmax * LRSF(i, j) / V);
                        
                        
                        A = syield / siiel;
                        AXY(i, j) = A;
                        % Count old yelding nodes
                        ynn = 0;
                        if(YNY0(i, j) > 0)
                            ynn = 1;
                            DSY(i, j) = SIIB(i, j) - syield;
                            ddd = ddd + DSY(i, j)^2;
                            ynpl = ynpl + 1;
                        end
                        % Update viscosity
                        if(A < 1)
                            % New viscosity for the basic node
                            etapl = dt * GGG(i, j) * A / (1 - A);
                            if(etapl < ETA0(i, j))
                                % Update plastic nodes
                                ETA5(i, j) = etapl^(1 - etawt) * ETA(i, j)^etawt;
                                YNY(i, j) = 1;
                                % Count yelding nodes
                                if(ynn == 0)
                                    DSY(i, j) = SIIB(i, j) - syield;
                                    ddd = ddd + DSY(i, j)^2;
                                    ynpl = ynpl + 1;
                                end
                            else
                                ETA5(i, j) = ETA0(i, j);
                            end
                        else
                            ETA5(i, j) = ETA0(i, j);
                        end
                    end
                end
            end
        end
        
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        % end measuring pressure cell time
        total_plasticity = total_plasticity + toc(t_plasticity);
        % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        % Compute Error
        DSYLSQ(iterstep) = 0;
        if(ynpl > 0);
            DSYLSQ(iterstep) = (ddd / ynpl)^0.5;
        end
        if(ynpl == 0)
            ETA = ETA0;
        end
        
        DSYERR = DSYLSQ(iterstep);
        
        % Adjust timestep
        dtpl = dt;
        % if (ynlast >= dtstep && ynpl > 0 && DSYLSQ(iterstep) > errmax && iterstep < niterglobal)
        %     dtpl = dt / dtkoef
        %     yn = 1;
        % end
        if (ynpl > 0 && iterstep < niterglobal && ynlast >= dtstep && (ynlast > ynlastmax || log10(DSYLSQ(iterstep) / DSYLSQ(iterstep - 1)) >= 0 || log10(DSYLSQ(iterstep) / DSYLSQ(iterstep - 1)) > log10(errmin / DSYLSQ(iterstep)) / (ynlastmax - ynlast)))
            dtpl = dt / dtkoef;
            yn = 1;
        end
        
        
        % Define displacement timesteps
        if(iterstep > 1)
            maxvxy0 = maxvxy;
        end
        maxvxs = max(max(abs(vxs)));
        maxvys = max(max(abs(vys)));
        maxvxy = ((max(max(vxs)) - min(min(vxs)))^2 + (max(max(vys)) - min(min(vys)))^2)^0.5;
        stpmaxcur = stpmax1;
        dtx = dt;
        if(dt > dx * stpmaxcur / maxvxs)
            dtx = dx / dtkoefv * stpmaxcur / maxvxs;
            yn = 1;
        end
        dty = dt;
        if(dt > dy * stpmaxcur / maxvys)
            dty = dy / dtkoefv * stpmaxcur / maxvys;
            yn = 1;
        end
        maxvxs = 0;
        maxvys = 0;
        for i = 1:1:Ny + 1
            for j = 1:1:Nx + 1
                if(yvx(i) >= upper_block && yvx(i) <= lower_block)
                    maxvxs = max(maxvxs, abs(vxs(i, j)));
                end
                if(yvy(i) >= upper_block && yvy(i) <= lower_block)
                    maxvys = max(maxvys, abs(vys(i, j)));
                end
            end
        end
        
        stpmaxcur = stpmax;
        if(dt > dx * stpmaxcur / maxvxs)
            dtx = dx / dtkoefv * stpmaxcur / maxvxs;
            yn = 1;
        end
        if(dt > dy * stpmaxcur / maxvys)
            dty = dy / dtkoefv * stpmaxcur / maxvys;
            yn = 1;
        end
        
        dtslip = 1e30;
        for i = 1:1:Ny
            for j = 1:1:Nx
                if(VSLIPB(i, j) > 0)
                    dtslip = min(dtslip, dx * stpmax / VSLIPB(i, j));
                end
            end
        end
        
        if(ynpl > 0 && dtslip < dt)
            yn = 1;
            dtslip = dtslip / dtkoefv;
        end
        
        
        % Chose minimal timestep
        if(yn > 0 && dt > dtmin)
            dtold = dt;
            dt = max(min([dtx dty dtpl dtslip dtlapusta]), dtmin);
            if(dt < dtold)
                ynlast = 0;
            end
        else
            yn = 0;
        end
        
        % Exit iterations
        ynstop = 0;
        % Velocity change ratio
        if(iterstep > 1)
            vratio = log10(maxvxy / maxvxy0);
        end
        if(yn == 0 && (ynpl == 0 || (DSYLSQ(iterstep) < errmin && iterstep > 1 && abs(vratio) < vratiomax)))
            ynstop = 1;
        else
            % Recomputing ETA
            for i = 1:1:Ny
                for j = 1:1:Nx
                    ETA(i, j) = max(min(ETA5(i, j), ETA0(i, j)), etamin);
                end
            end
            % Save current viscosity
            ETA50 = ETA;
            YNY0 = YNY;
            OM = OM5;
        end
        
        
        % Exit iteration
        if(ynstop == 1)
            break;
        end
        ynlast = ynlast + 1;
    end
    
    % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    % end measuring marker time
    total_time_global = total_time_global + toc(t_global_iterations);
    % """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    % /////////////////////////////////////////////////////////////////////////////////////// 
    % end of loop through globel iterations
    % /////////////////////////////////////////////////////////////////////////////////////// 
    
    % Mark dt decrease
    if(dt00 > dt)
        yndtdecrease = 1;
    end
    
    % Save current viscosity
    ETA00 = ETA50;
    OM0 = OM;
    % Recheck displacement timestep
    dt = min([dtx dty dt dtlapusta]);
    
    % Compute strain rate,  stress and stress change
    ESP = zeros(Ny, Nx);
    EXY = zeros(Ny, Nx);
    SXY = zeros(Ny, Nx);
    DSXY = zeros(Ny, Nx);
    EXX = zeros(Ny1, Nx1);
    SXX = zeros(Ny1, Nx1);
    DSXX = zeros(Ny1, Nx1);
    EYY = zeros(Ny1, Nx1);
    SYY = zeros(Ny1, Nx1);
    DSYY = zeros(Ny1, Nx1);
    EII = zeros(Ny1, Nx1);
    EIIVP = zeros(Ny1, Nx1);
    SII = zeros(Ny1, Nx1);
    DSII = zeros(Ny1, Nx1);
    DIS = zeros(Ny1, Nx1);
    
    % /////////////////////////////////////////////////////////////////////////////////////// 
    % third layer loop (timestep -> Ny -> Nx), 
    % /////////////////////////////////////////////////////////////////////////////////////// 
    
    % Process internal basic nodes
    for i = 1:1:Ny
        for j = 1:1:Nx
            % ESP = 1 / 2(dVy / dx - dVx / dy),  EXY,  SXY,  DSXY
            ESP(i, j) = 1 / 2 * ((vys(i, j + 1) - vys(i, j)) / dx - (vxs(i + 1, j) - vxs(i, j)) / dy);
            EXY(i, j) = 1 / 2 * ((vxs(i + 1, j) - vxs(i, j)) / dy + (vys(i, j + 1) - vys(i, j)) / dx);
            KXY = dt * GGG(i, j) / (dt * GGG(i, j) + ETA(i, j));
            SXY(i, j) = 2 * ETA(i, j) * EXY(i, j) * KXY + SXY0(i, j) * (1 - KXY);
            DSXY(i, j) = SXY(i, j) - SXY0(i, j);
        end
    end
    
    % /////////////////////////////////////////////////////////////////////////////////////// 
    % very similar loop as above,  maybe merge
    % third layer loop (timestep -> Ny -> Nx), 
    % /////////////////////////////////////////////////////////////////////////////////////// 
    
    % Process pressure cells
    for i = 2:1:Ny
        for j = 2:1:Nx
            % EXX,  SXX,  DSXX
            EXX(i, j) = (2 * (vxs(i, j) - vxs(i, j - 1)) / dx - (vys(i, j) - vys(i - 1, j)) / dy) / 3;
            EYY(i, j) = (2 * (vys(i, j) - vys(i - 1, j)) / dy - (vxs(i, j) - vxs(i, j - 1)) / dx) / 3;
            KXX = dt * GGGP(i, j) / (dt * GGGP(i, j) + ETAP(i, j));
            SXX(i, j) = 2 * ETAP(i, j) * EXX(i, j) * KXX + SXX0(i, j) * (1 - KXX);
            SYY(i, j) = 2 * ETAP(i, j) * EYY(i, j) * KXX + SYY0(i, j) * (1 - KXX);
            DSXX(i, j) = SXX(i, j) - SXX0(i, j);
            DSYY(i, j) = SYY(i, j) - SYY0(i, j);
        end
    end
    % Compute stress and strain rate invariants
    % Process pressure cells
    for i = 2:1:Ny
        for j = 2:1:Nx
            % EXY term is averaged from four surrounding basic nodes
            EXY1 = (EXY(i, j) + EXY(i - 1, j) + EXY(i, j - 1) + EXY(i - 1, j - 1)) / 4;
            EII(i, j) = (0.5 * (EXX(i, j)^2 + EYY(i, j)^2) + EXY1^2)^0.5;
            % Second strain rate invariant SII
            % SXY term is averaged from four surrounding basic nodes
            SXY1 = (SXY(i, j) + SXY(i - 1, j) + SXY(i, j - 1) + SXY(i - 1, j - 1)) / 4;
            SII(i, j) = (0.5 * (SXX(i, j)^2 + SYY(i, j)^2) + SXY1^2)^0.5;
        end
    end
    
    % Compute stress and strain rate invariants and dissipation
    % Process pressure cells
    for i = 2:1:Ny
        for j = 2:1:Nx
            % EXY term is averaged from four surrounding basic nodes
            EXY2 = (EXY(i, j)^2 + EXY(i - 1, j)^2 + EXY(i, j - 1)^2 + EXY(i - 1, j - 1)^2) / 4;
            EII(i, j) = (EXX(i, j)^2 + EXY2)^0.5;
            EXYVP2 = ((SXY(i, j) / 2 / ETA(i, j))^2 + (SXY(i - 1, j) / 2 / ETA(i - 1, j))^2 + (SXY(i, j - 1) / 2 / ETA(i, j - 1))^2 + (SXY(i - 1, j - 1) / 2 / ETA(i - 1, j - 1))^2) / 4;
            EIIVP(i, j) = (0.5 * ((SXX(i, j) / 2 / ETAP(i, j))^2 + (SYY(i, j) / 2 / ETAP(i, j))^2) + EXYVP2)^0.5;
            % Second strain rate invariant SII
            % SXY term is averaged from four surrounding basic nodes
            SXY2 = (SXY(i, j)^2 + SXY(i - 1, j)^2 + SXY(i, j - 1)^2 + SXY(i - 1, j - 1)^2) / 4;
            SII(i, j) = (0.5 * (SXX(i, j)^2 + SYY(i, j)^2) + SXY2)^0.5;
            
            % Dissipation
            DISXY = (SXY(i, j)^2 / 2 / ETA(i, j) + SXY(i - 1, j)^2 / 2 / ETA(i - 1, j) + SXY(i, j - 1)^2 / 2 / ETA(i, j - 1) + SXY(i - 1, j - 1)^2 / 2 / ETA(i - 1, j - 1)) / 4;
            DIS(i, j) = SXX(i, j)^2 / 2 / ETAP(i, j) + SYY(i, j)^2 / 2 / ETAP(i, j) + 2 * DISXY;
            
            
            if(i < Ny)
                pt_ave(i, j)  = (pt(i, j) + pt(i + 1, j)) / 2;
                pf_ave(i, j)  = (pf(i, j) + pf(i + 1, j)) / 2;
                PT0_ave      = (PT0(i, j) + PT0(i + 1, j)) / 2;
                PF0_ave      = (PF0(i, j) + PF0(i + 1, j)) / 2;
            else
                pt_ave(i, j)  = pt(i, j);
                pf_ave(i, j)  = pf(i, j);
                PT0_ave      = PT0(i, j);
                PF0_ave      = PF0(i, j);
            end
            
            % Compute elastic and viscous compaction
            VIS_COMP(i, j) = (pt_ave(i, j) - pf_ave(i, j)) / (ETAB(i, j) * (1 - POR(i, j)));
            % Drained compressibility
            BETTADRAINED = (1 / GGGB(i, j) + BETTASOLID) / (1 - POR(i, j));
            % Biott - Willis koefficient
            KBW = 1 - BETTASOLID / BETTADRAINED;
            EL_DECOM(i, j) = BETTADRAINED * (pt_ave(i, j) - PT0_ave-KBW * pf_ave(i, j) + KBW * PF0_ave) / dt;
            
        end
    end
    
    
    % Runge-Kutta velocity,  spin array
    vxm = zeros(4, 1);
    vym = zeros(4, 1);
    spm = zeros(4, 1);
    % Move markers by nodal velocity field
    for m = 1:1:marknum
        % Save marker position
        xold = xm(m);
        yold = ym(m);
        for rk = 1:1:4
            % vx - velocity interpolation
            % [i, j] -------- [i, j + 1]
            %   |                |
            %   |    o m         |
            %   |                |
            % [i + 1, j] ------- [i + 1, j + 1]
            % Indexes and distances
            j = fix(xm(m) / dx) + 1;
            i = fix((ym(m) + dy / 2) / dy) + 1;
            if(j < 1)
                j = 1;
            elseif (j > Nx - 1)
                j = Nx - 1;
            end
            if(i < 1)
                i = 1;
            elseif (i > Ny)
                i = Ny;
            end
            %Distances
            dxm = (xm(m) - xvx(j)) / dx;
            dym = (ym(m) - yvx(i)) / dy;
            % Weights
            wtmij = (1 - dxm) * (1 - dym);
            wtmi1j = (1 - dxm) * (dym);
            wtmij1 = (dxm) * (1 - dym);
            wtmi1j1 = (dxm) * (dym);
            % Interpolation
            vxm(rk) = vxs(i, j) * wtmij + vxs(i + 1, j) * wtmi1j + vxs(i, j + 1) * wtmij1 + vxs(i + 1, j + 1) * wtmi1j1;
            
            % vy - velocity interpolation
            % [i, j] -------- [i, j + 1]
            %   |                |
            %   |    o m         |
            %   |                |
            % [i + 1, j] ------- [i + 1, j + 1]
            % Indexes and distances
            j = fix((xm(m) + dx / 2) / dx) + 1;
            i = fix(ym(m) / dy) + 1;
            if(j < 1)
                j = 1;
            elseif (j > Nx)
                j = Nx;
            end
            if(i < 1)
                i = 1;
            elseif (i > Ny - 1)
                i = Ny - 1;
            end
            %Distances
            dxm = (xm(m) - xvy(j)) / dx;
            dym = (ym(m) - yvy(i)) / dy;
            % Weights
            wtmij = (1 - dxm) * (1 - dym);
            wtmi1j = (1 - dxm) * (dym);
            wtmij1 = (dxm) * (1 - dym);
            wtmi1j1 = (dxm) * (dym);
            % Interpolation
            vym(rk) = vys(i, j) * wtmij + vys(i + 1, j) * wtmi1j + vys(i, j + 1) * wtmij1 + vys(i + 1, j + 1) * wtmi1j1;
            
            
            % ESP = 1 / 2(dVy / dx - dVx / dy) interpolation
            % [i, j] -------- [i, j + 1]
            %   |                |
            %   |    o m         |
            %   |                |
            % [i + 1, j] ------- [i + 1, j + 1]
            % Indexes and distances
            j = fix((xm(m)) / dx) + 1;
            i = fix((ym(m)) / dy) + 1;
            if(j < 1)
                j = 1;
            elseif (j > Nx - 1)
                j = Nx - 1;
            end
            if(i < 1)
                i = 1;
            elseif (i > Ny - 1)
                i = Ny - 1;
            end
            
            %Distances
            dxm = (xm(m) - x(j)) / dx;
            dym = (ym(m) - y(i)) / dy;
            % Weights
            wtmij = (1 - dxm) * (1 - dym);
            wtmi1j = (1 - dxm) * (dym);
            wtmij1 = (dxm) * (1 - dym);
            wtmi1j1 = (dxm) * (dym);
            % Interpolation ESP = 1 / 2(dVy / dx - dVx / dy) for the marker
            spm(rk) = ESP(i, j) * wtmij + ESP(i + 1, j) * wtmi1j + ESP(i, j + 1) * wtmij1 + ESP(i + 1, j + 1) * wtmi1j1;
            
            
            % Moving between A, B, C, D points
            if(rk < 3)
                % Moving A -  > B and A -  > C
                xm(m) = xold + vxm(rk) * dt / 2;
                ym(m) = yold + vym(rk) * dt / 2;
            elseif(rk == 3)
                % Moving A -  > D
                xm(m) = xold + vxm(rk) * dt;
                ym(m) = yold + vym(rk) * dt;
            end
        end
        % Compute effective velocity,  rotation rate
        vxeff = (vxm(1) + 2 * vxm(2) + 2 * vxm(3) + vxm(4)) / 6;
        vyeff = (vym(1) + 2 * vym(2) + 2 * vym(3) + vym(4)) / 6;
        speff = spm(1);
        
        % Rotate stress on marker according to its spin
        % Compute amount of rotation from spin rate:
        % Espin = 1 / 2(dvy / dx - dvx / dy) i.e. positive for clockwise rotation
        % (when x axis is directed rightward and y axis is directed downward)
        dspeff = speff * dt;
        % Save old stresses
        msxxold = sxxm(m);
        msyyold = syym(m);
        msxyold = sxym(m);
        sxym(m) = 0.5 * (msxxold - msyyold) * sin(2 * dspeff) + msxyold * cos(2 * dspeff);
        sxxm(m) = msxxold * (cos(dspeff))^2 + msyyold * (sin(dspeff))^2 - msxyold * sin(2 * dspeff);
        syym(m) = msxxold * (sin(dspeff))^2 + msyyold * (cos(dspeff))^2 + msxyold * sin(2 * dspeff);
        
        % Move markers
        xm(m) = xold + vxeff * dt;
        ym(m) = yold + vyeff * dt;
        
        % Recycling
        if(xm(m) < 0)
            xm(m) = xm(m) + xsize;
        end
        if(xm(m) > xsize)
            xm(m) = xm(m) - xsize;
        end
    end
    
    % Update timesum
    timesum = timesum + dt;
    timesumcur(timestep) = timesum;
    timemyr = timesum / (1e+6 * 365.25 * 24 * 3600);
    timeyr = timesum / (365.25 * 24 * 3600);
    timehr = timesum / (3600);
    dtcur(timestep) = dt;
    
    maxvxsmod(timestep) = -1e+30;
    minvxsmod(timestep) = 1e+30;
    maxvysmod(timestep) = -1e+30;
    minvysmod(timestep) = 1e+30;
    
    % /////////////////////////////////////////////////////////////////////////////////////// 
    % third layer loop (timestep -> Ny -> Nx), 
    % /////////////////////////////////////////////////////////////////////////////////////// 
    
    % Vx
    for i = 2:1:Ny
        for j = 1:1:Nx
            if(RHOX(i, j) > 2000)
                maxvxsmod(timestep) = max(maxvxsmod(timestep), vxs(i, j));
                minvxsmod(timestep) = min(minvxsmod(timestep), vxs(i, j));
            end
        end
    end
    % Vy
    for i = 1:1:Ny
        for j = 2:1:Nx
            if(RHOY(i, j) > 2000)
                maxvysmod(timestep) = max(maxvysmod(timestep), vys(i, j));
                minvysmod(timestep) = min(minvysmod(timestep), vys(i, j));
            end
        end
    end
    
    % Update VX0
    VX0 = vxs;
    for i = 1:1:Ny1
        for j = 1:1:Nx1
            if(PORX(i, j) > 0)
                VXF0(i, j) = vxs(i, j) + vxD(i, j) / PORX(i, j);
            end
        end
    end
    % Update VY0
    VY0 = vys;
    for i = 1:1:Ny1
        for j = 1:1:Nx1
            if(PORY(i, j) > 0)
                VYF0(i, j) = vys(i, j) + vyD(i, j) / PORY(i, j);
            end
        end
    end
    % Update SXX0
    SXX0 = SXX;
    % Update SYY0
    SYY0 = SYY;
    % Update SXY0
    SXY0 = SXY;
    % Update PTF0
    PTF0 = (pt - pf);
    PT0 = pt;
    PF0 = pf;
    
    % /////////////////////////////////////////////////////////////////////////////////////// 
    % output
    % /////////////////////////////////////////////////////////////////////////////////////// 
    
    fprintf(' ==================================== \n');
    fprintf('total time:        %.12E sec \n', timesum);
    fprintf('time step:         %.12E sec \n', dt);
    fprintf('Vslip max:         %.12E m / s \n', Vmax);
    fprintf('iter - iterations:   %d \n', iterstep);
    fprintf('global - iterations: %d \n', ynlast);
    
    if seismic_cycles_data
        
        if (timesum == dt)
            data_save = x;
            fileID    = fopen('x_fault.txt', 'a');
            TP_write  = fprintf (fileID, '%.3E    \n', data_save);
            fclose(fileID);
        end
        
        if timesum_plus  <  timesum
            
            % ==========  save slip rate
            data_save = [timesum dt VSLIPB(line_fault, :)];
            fid = fopen('EVO_Vslip.txt', 'a');
            fprintf(fid, '%.6E ', data_save);
            fprintf(fid, '\n');
            fclose(fid);
            clear data_save
            
            % ==========  save viscosity
            data_save = [timesum dt ETA(line_fault, :)];
            fid = fopen('EVO_viscosity.txt', 'a');
            fprintf(fid, '%.6E ', data_save);
            fprintf(fid, '\n');
            fclose(fid);
            clear data_save
            
            % ==========  save fluid pressure
            data_save = [timesum dt pf_ave(line_fault, :)];
            fid = fopen('EVO_press_flu.txt', 'a');
            fprintf(fid, '%.6E ', data_save);
            fprintf(fid, '\n');
            fclose(fid);
            clear data_save
            
            % ==========  save effective pressure
            P_diff = pt_ave-pf_ave;
            data_save = [timesum dt P_diff(line_fault, :)];
            fid = fopen('EVO_press_eff.txt', 'a');
            fprintf(fid, '%.6E ', data_save);
            fprintf(fid, '\n');
            fclose(fid);
            clear data_save
            
            % ==========  save SigmaY
            data_save = [timesum dt SigmaY(line_fault, :)];
            fid = fopen('EVO_SigmaY.txt', 'a');
            fprintf(fid, '%.9E ', data_save);
            fprintf(fid, '\n');
            fclose(fid);
            clear data_save
            
            % ==========  save SII
            data_save = [timesum dt SII_fault(line_fault, :)];
            fid = fopen('EVO_Sii.txt', 'a');
            fprintf(fid, '%.9E ', data_save);
            fprintf(fid, '\n');
            fclose(fid);
            clear data_save
            
            % ==========  save Theta
            data_save = [timesum dt OM(line_fault, :)];
            fid = fopen('EVO_Theta.txt', 'a');
            fprintf(fid, '%.9E ', data_save);
            fprintf(fid, '\n');
            fclose(fid);
            clear data_save
            
            % ==========  save viscous compaction
            data_save = [timesum dt VIS_COMP(line_fault, :)];
            fid = fopen('EVO_Visc_comp.txt', 'a');
            fprintf(fid, '%.9E ', data_save);
            fprintf(fid, '\n');
            fclose(fid);
            clear data_save
            
            % ==========  save elastic compaction
            data_save = [timesum dt EL_DECOM(line_fault, :)];
            fid = fopen('EVO_Elast_comp.txt', 'a');
            fprintf(fid, '%.9E ', data_save);
            fprintf(fid, '\n');
            fclose(fid);
            clear data_save
            
            % ==========  save time,  dt,  vmax
            Vmax = max(max(VSLIPB));
            fileID    = fopen('EVO_data.txt', 'a');
            TP_write  = fprintf (fileID, '%.12E  %.12E  %.12E  %d    %d \n', timesum, dt, Vmax, ynlast, iterstep);
            fclose(fileID);
            
            timesum_plus = timesum;
        end
    end
    
    if(fix(timestep / savematstep) * savematstep == timestep)
        namemat    =  [nname, num2str(timestep)];
        save(namemat);
        fdata = fopen('file.txt', 'wt');
        fprintf(fdata, '%d', timestep);
        fclose(fdata);
    end
    
end

% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
% stop measuring initialisation time
iteration_time = toc(t_iteration);
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

% /////////////////////////////////////////////////////////////////////////////////////// 
% end of loop through timesteps
% /////////////////////////////////////////////////////////////////////////////////////// 

% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
% end measuring time
time_elapsed = toc(t_start);
time_output = fopen('time.txt', 'a');

fprintf(time_output,  'total time: \n');
fprintf(time_output,  '%d seconds',  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Initialisation time: \n');
fprintf(time_output,  '%d seconds',  initialisation_time);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * initialisation_time  /  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Total iteration time: \n');
fprintf(time_output,  '%d seconds',  iteration_time);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * iteration_time  /  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Average iteration time: \n');
fprintf(time_output,  '%d seconds',  iteration_time  /  num_timesteps);
fprintf(time_output,  '\n -------------------------- \n');

fprintf(time_output,  'Total markers time: \n');
fprintf(time_output,  '%d seconds',  total_time_markers);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * total_time_markers  /  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Total nodes time: \n');
fprintf(time_output,  '%d seconds',  total_nodes);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * total_nodes  /  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Total global iteration time: \n');
fprintf(time_output,  '%d seconds',  total_time_global);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * total_time_global  /  time_elapsed);
fprintf(time_output,  '\n -------------------------- \n');

fprintf(time_output,  'Total plastic strain time: \n');
fprintf(time_output,  '%d seconds',  total_plastic_strain);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * total_plastic_strain  /  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Total viscoplastic time: \n');
fprintf(time_output,  '%d seconds',  total_viscoplastic);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * total_viscoplastic  /  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Total global matrices time: \n');
fprintf(time_output,  '%d seconds',  total_global_matrices);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * total_global_matrices  /  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Total solve matrices time: \n');
fprintf(time_output,  '%d seconds',  total_solve);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * total_solve  /  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Total reload solution time: \n');
fprintf(time_output,  '%d seconds',  total_reload_solution);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * total_reload_solution  /  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Total internal nodes time: \n');
fprintf(time_output,  '%d seconds',  total_internal_nodes);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * total_internal_nodes  /  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Total pressure cell time: \n');
fprintf(time_output,  '%d seconds',  total_pressure_cell);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * total_pressure_cell  /  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Total pressure cell2 time: \n');
fprintf(time_output,  '%d seconds',  total_pressure_cell2);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * total_pressure_cell2  /  time_elapsed);
fprintf(time_output,  '\n\n');

fprintf(time_output,  'Total plasticity time: \n');
fprintf(time_output,  '%d seconds',  total_plasticity);
fprintf(time_output,  '\n');
fprintf(time_output,  '%.4f %%',  100 * total_plasticity  /  time_elapsed);
fprintf(time_output,  '\n\n');

fclose(time_output);
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""