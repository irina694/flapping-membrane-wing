function[g] = membrane_flapping_fun(x,plot_decide)

g_fail = [1,1,1];  % return these values for the objective functions if something fails
N_elements_fail = 6000;  % if the grid has more than this # of elements, don't bother to evaluate design, g = g_fail

%% Wing geometry:

wing.root_chord = 0.16; % root chord of wing (m)
wing.length = .4; % length of wing (m)
wing.x_offset = -.04; % x-center of rotation (m)
wing.y_offset = 0; % y-center of rotation (m)
wing.camber = .02; % parabolic camber as a fraction of local chord
wing.twist = 0; % built-in twist at wing tip (degrees)

wing.y_stations = [0;.25;.5;.75;.9;.95;1]*wing.length; % spanwise location of knots used to define wing planform
wing.LE_knots = [0;9.526E-4;4.0192E-3;.0102;.0169;.0206;.03]; % knots used to define LE
wing.TE_knots = [wing.root_chord;.1571;.1479;.1296;.1093;.0981;.07]; % knots used to define TE

loads.tip_loss = 0.7; % governs losses due to tip effects
loads.cd_o = 0.05; loads.cd_pi2 = 2; % viscous drag coefficients at AOA = 0, pi/2

%% Kinematics:

kin.omega = 40; % flapping frequency (rad/s)

kin.z_m = 0; % amplitude of sinusoidal z-plunging motion (m)
kin.z_shift = 0; % phase shift of sinusoidal z-plunging motion (rad)
kin.z_DC = 0; % offset of sinusoidal z-plunging motion (m)

kin.f_m = pi/6; % amplitude of sinusoidal flapping motion (rad)
kin.f_shift = 0; % phase shift of sinusoidal flapping motion (rad)
kin.f_DC = 0; % offset of sinusoidal flapping motion (rad)

kin.th_m = 0; % amplitude of sinusoidal pitching motion (rad)
kin.th_shift = 0; % phase shift of sinusoidal pitching motion (rad)
kin.th_DC = pi/45; % offset of sinusoidal pitching motion (rad)

kin.periodic = true; % tells the solver whether or not the same prescribed kinematics may be expected during each cycle

%% time step information:

time.N_per = 100; % how many time steps per flapping cycle?
time.N_cycles = 5; % how many flapping cycles (want to time march until initial transients die out)?
time.skip = 10; % only plot output every time.skip time steps

time.t_step = 2*pi/kin.omega/time.N_per;
time.N_steps = time.N_per*time.N_cycles;
time.t = (0:time.t_step:time.t_step*time.N_steps)';

%% Free-stream conditions:

setup.U = 10*ones(length(time.t),1); % free-stream velocity at each time step (m/s)
setup.Up = 0*ones(length(time.t),1); % time-derivative of free-stream velocity at each time step (m/s/s)
setup.rho = 1.225; % flow density (kg/m^3)

%% Wing topology

[connectivity,connectivity.error] = create_cell(x); % create connectivity
connectivity.gage = .02; % seeding of control points for mesh generation (m)

if connectivity.error == 0
    
    %% Wing structure
    
    mat.membrane.E = 2E6; % membrane skin modulus (Pa)
    mat.membrane.poisson = 0.5; % membrane skin Poisson ratio
    mat.membrane.thickness = .1/1000; % membrane skin thickness (m)
    mat.membrane.pre_stress = 10; % membrane skin pre-stress (N/m)
    mat.membrane.rho = 1200; % membrane skin density (kg/m^3)
    
    mat.LE.thickness = 2/1000; % thickness of the leading edge battens (m)
    mat.LE.base = 5/1000; % width of the leading edge battens (m)
    mat.LE.E = 300E9; % modulus of the leading edge battens (Pa)
    mat.LE.rho = 1600; % density of the leading edge battens (kg/m^3)
    
    % the trailing edge beams are removed below, so these numbers don't matter
    mat.TE.thickness = .2/1000; % thickness of the trailing edge battens (m) 
    mat.TE.base = 3/1000; % width of the trailing edge battens (m) 
    mat.TE.E = 300E9; % modulus of the trailing edge battens (Pa) 
    mat.TE.rho = 1600; % density of the trailing edge battens (kg/m^3)
    
    mat.TIP.thickness = .8/1000; % thickness of the wingtip battens (m) 
    mat.TIP.base = 3/1000; % width of the wingtip battens (m)  
    mat.TIP.E = 300E9; % modulus of the wingtip battens (Pa)  
    mat.TIP.rho = 1600; % density of the wingtip battens (kg/m^3)
    
    mat.ROOT.thickness = .8/1000; % thickness of the root battens (m)  
    mat.ROOT.base = 3/1000; % width of the root battens (m)   
    mat.ROOT.E = 300E9; % modulus of the root battens (Pa) 
    mat.ROOT.rho = 1600; % density of the root battens (kg/m^3)
    
    mat.INT.thickness = .8/1000; % thickness of the internal battens (m)  
    mat.INT.base = 3/1000; % width of the internal battens (m)  
    mat.INT.E = 300E9; % modulus of the internal battens (Pa) 
    mat.INT.rho = 1600; % density of the internal battens (kg/m^3)
    
    mat.damping = 10; % mass-proportional damping of the entire system
    
    FEA.N_modes = 20; % number of finite element modes
    
    %% Alpha integration parameters
    
    solve.alpha = -.25; solve.beta = 0.75;
    
    %% Number of Glauert expansions:
    
    aero.N_h = 10;
    
    %% Number of inflow expansions:
    
    inflow.N = 6;
    
    %% Gauss points:
    
    gauss.N_gauss = 20; % number of guass points at each span station
    gauss.N_stations = 20; % number of span stations
    
    %% Kinematic motions:
    
    kin.z = kin.z_m*sin(kin.omega*time.t + kin.z_shift) + kin.z_DC; % kinematic motions at each time step
    kin.zp = kin.omega*kin.z_m*cos(kin.omega*time.t + kin.z_shift);
    kin.zpp = -kin.omega*kin.omega*kin.z_m*sin(kin.omega*time.t + kin.z_shift);
    
    kin.f = kin.f_m*sin(kin.omega*time.t + kin.f_shift) + kin.f_DC;
    kin.fp = kin.omega*kin.f_m*cos(kin.omega*time.t + kin.f_shift);
    kin.fpp = -kin.omega*kin.omega*kin.f_m*sin(kin.omega*time.t + kin.f_shift);
    
    kin.th = kin.th_m*sin(kin.omega*time.t + kin.th_shift) + kin.th_DC;
    kin.thp = kin.omega*kin.th_m*cos(kin.omega*time.t + kin.th_shift);
    kin.thpp = -kin.omega*kin.omega*kin.th_m*sin(kin.omega*time.t + kin.th_shift);
    
    %% Define wing geometry:
    
    [wing.x,wing.y,wing.z,wing.LE_pp,wing.TE_pp,gauss.weight,gauss.theta,gauss.x,gauss.y,gauss.zo,gauss.LE,gauss.b,...
        connectivity,FEA.seed_node,FEA.seed_edge,FEA.x,FEA.y,FEA.z,FEA.trielements,FEA.beam_elements,FEA.beam_ID] = ...
        wing_geom2(wing.root_chord,wing.length,wing.camber,wing.twist,wing.y_stations,wing.LE_knots,wing.TE_knots,gauss.N_gauss,gauss.N_stations,connectivity);
    FEA.beam_elements(FEA.beam_ID == 'S',:) = []; % remove the beam elements at the trailing edge
    FEA.beam_ID(FEA.beam_ID == 'S') = [];
    
    if length(FEA.trielements) < N_elements_fail % only test the system if it isn't too dense
        
        loads.correct_station = 1-loads.tip_loss*exp((2*2*wing.length/wing.root_chord)*(gauss.y/wing.length-1)); % load correction factors at each span station due to tip losses
        loads.correct_total = sparse(1:gauss.N_gauss*gauss.N_stations,1:gauss.N_gauss*gauss.N_stations,reshape(repmat(loads.correct_station,gauss.N_gauss,1),gauss.N_gauss*gauss.N_stations,1));
        
        %% Glauert expansion and contraction matrices:
        
        interp.T_hn_h = repmat(gauss.weight',aero.N_h+1,1).*cos(repmat([0:aero.N_h]',1,gauss.N_gauss).*repmat(gauss.theta',aero.N_h+1,1));
        interp.T_hn_h(1,:) = interp.T_hn_h(1,:)/pi; interp.T_hn_h(2:end,:) = interp.T_hn_h(2:end,:)*2/pi;
        interp.T_hn_h = kron(eye(gauss.N_stations),interp.T_hn_h);
        
        interp.P_map = zeros(gauss.N_gauss,aero.N_h+1);
        interp.P_map(:,1) = (1-cos(gauss.theta))./sin(gauss.theta);
        interp.P_map(:,2:aero.N_h+1) = sin(repmat(1:aero.N_h,gauss.N_gauss,1).*repmat(gauss.theta,1,aero.N_h));
        interp.P_map = kron(eye(gauss.N_stations),interp.P_map);
        
        interp.h_map = cos(repmat(0:aero.N_h,gauss.N_gauss,1).*repmat(gauss.theta,1,aero.N_h+1));
        interp.h_map = kron(eye(gauss.N_stations),interp.h_map);
        
        interp.hx_map = zeros(gauss.N_gauss,aero.N_h+1);
        interp.hx_map(:,2:aero.N_h+1) = repmat(1:aero.N_h,gauss.N_gauss,1).*sin(repmat(1:aero.N_h,gauss.N_gauss,1).*repmat(gauss.theta,1,aero.N_h))./sin(repmat(gauss.theta,1,aero.N_h));
        interp.hx_map = kron(diag(1./gauss.b),interp.hx_map);
        
        %% Aero loading terms:
        
        aero.uo = repmat(setup.U.*cos(kin.th),1,gauss.N_stations) - (repmat(gauss.y,time.N_steps+1,1) + wing.y_offset).*repmat(kin.thp.*sin(kin.f),1,gauss.N_stations) + repmat(kin.zp.*sin(kin.th),1,gauss.N_stations);
        aero.vo = -repmat(gauss.y,time.N_steps+1,1).*repmat(kin.fp,1,gauss.N_stations) + repmat((setup.U.*sin(kin.f+kin.th))/2 - kin.fp*wing.y_offset - (kin.zp.*cos(kin.f+kin.th))/2  - (kin.zp.*cos(kin.f-kin.th))/2 - (setup.U.*sin(kin.f-kin.th))/2,1,gauss.N_stations) + repmat(kin.thp.*cos(kin.f),1,gauss.N_stations).*(wing.x_offset + repmat(gauss.LE+gauss.b,time.N_steps+1,1));
        aero.v1 = repmat(gauss.b,time.N_steps+1,1).*repmat(kin.thp.*cos(kin.f),1,gauss.N_stations);
        aero.uo_dot = repmat(setup.Up.*cos(kin.th),1,gauss.N_stations) - repmat(setup.U.*sin(kin.th).*kin.thp,1,gauss.N_stations) - (repmat(gauss.y,time.N_steps+1,1) + wing.y_offset).*repmat(kin.thpp.*sin(kin.f),1,gauss.N_stations) - (repmat(gauss.y,time.N_steps+1,1) + wing.y_offset).*repmat(kin.thp.*cos(kin.f).*kin.fp,1,gauss.N_stations) + repmat(kin.zpp.*sin(kin.th),1,gauss.N_stations) + repmat(kin.zp.*cos(kin.th).*kin.thp,1,gauss.N_stations);
        aero.vo_dot = -repmat(gauss.y,time.N_steps+1,1).*repmat(kin.fpp,1,gauss.N_stations) + repmat((setup.Up.*sin(kin.f+kin.th))/2 + (setup.U.*cos(kin.f+kin.th).*(kin.fp+kin.thp))/2 - kin.fpp*wing.y_offset - (kin.zpp.*cos(kin.f+kin.th))/2 + (kin.zp.*sin(kin.f+kin.th).*(kin.fp+kin.thp))/2 - (kin.zpp.*cos(kin.f-kin.th))/2 + (kin.zp.*sin(kin.f-kin.th).*(kin.fp-kin.thp))/2 - (setup.Up.*sin(kin.f-kin.th))/2 - (setup.U.*cos(kin.f-kin.th).*(kin.fp-kin.thp))/2,1,gauss.N_stations) + repmat(kin.thpp.*cos(kin.f) - kin.thp.*sin(kin.f).*kin.fp,1,gauss.N_stations).*(wing.x_offset + repmat(gauss.LE+gauss.b,time.N_steps+1,1));
        aero.v1_dot = repmat(gauss.b,time.N_steps+1,1).*repmat(kin.thpp.*cos(kin.f)-kin.thp.*sin(kin.f).*kin.fp,1,gauss.N_stations);
        
        aero.vn = zeros((aero.N_h+1)*gauss.N_stations,length(time.t)); aero.vn(1:(aero.N_h+1):size(aero.vn,1),:) = aero.vo'; aero.vn(2:(aero.N_h+1):size(aero.vn,1),:) = aero.v1';
        aero.vn_dot = zeros((aero.N_h+1)*gauss.N_stations,length(time.t)); aero.vn_dot(1:(aero.N_h+1):size(aero.vn_dot,1),:) = aero.vo_dot'; aero.vn_dot(2:(aero.N_h+1):size(aero.vn_dot,1),:) = aero.v1_dot';
        
        aero.K = zeros(aero.N_h+1,aero.N_h+1);
        if mod(aero.N_h,2) == 0
            aero.K(1,[2:2:aero.N_h+1]) = [1:2:aero.N_h];
            for i = 1:aero.N_h/2
                aero.K = aero.K + diag([0,i*4:2:aero.N_h*2],(i-1)*2+1);
            end
        else
            aero.K(1,[2:2:aero.N_h+1]) = [1:2:aero.N_h+1];
            for i = 1:aero.N_h/2
                aero.K = aero.K + diag([0,i*4:2:aero.N_h*2],(i-1)*2+1);
            end
        end
        aero.K = kron(diag(1./gauss.b),aero.K);
        
        aero.W = zeros(aero.N_h+1,aero.N_h+1); aero.W(2,1) = 2; aero.W(2,3) = -1;
        aero.W = aero.W + diag([0,1./[2:1:aero.N_h]],-1);
        aero.W = aero.W - diag([0,0,1./[2:1:aero.N_h-1]],1);
        aero.W = kron(diag(gauss.b/2),aero.W);
        
        %% Inflow matrices:
        
        inflow.b = (((-1).^((1:inflow.N-1)-1)).*factorial(inflow.N+(1:inflow.N-1)-1)./factorial(inflow.N-(1:inflow.N-1)-1)./(factorial(1:inflow.N-1).^2))';
        inflow.b(end+1,1) = (-1)^(inflow.N+1);
        inflow.c = 2./(1:inflow.N)';
        inflow.d = 0*inflow.c; inflow.d(1) = 0.5;
        inflow.D = diag(1./(2:inflow.N)/2,-1)-diag(1./(1:inflow.N-1)/2,1);
        inflow.A = inflow.D + inflow.d*inflow.b' + inflow.c*inflow.d' + .5*inflow.c*inflow.b';
        
        inflow.A = kron(speye(gauss.N_stations),inflow.A);
        inflow.Q = kron(speye(gauss.N_stations),[1;zeros(aero.N_h,1)])*kron(speye(gauss.N_stations),.5*inflow.b');
        inflow.Rw = kron(speye(gauss.N_stations),inflow.c)*kron(speye(gauss.N_stations),[1;.5;zeros(aero.N_h-1,1)])';
        
        %% FEA parameters:
        
        FEA.N_nodes = length(FEA.x);
        FEA.GlobalDOF = FEA.N_nodes*6;
        FEA.DOF = reshape(1:FEA.GlobalDOF,6,FEA.N_nodes)';
        
        for i = 1:length(FEA.trielements)
            FEA.membrane_asmb(i,:) = [FEA.DOF(FEA.trielements(i,1),:),FEA.DOF(FEA.trielements(i,2),:),FEA.DOF(FEA.trielements(i,3),:)];
        end
        for i = 1:length(FEA.beam_elements)
            FEA.beam_asmb(i,:) = [FEA.DOF(FEA.beam_elements(i,1),:),FEA.DOF(FEA.beam_elements(i,2),:)];
        end
        
        %% FEA boundary conditions
        
        FEA.bc = FEA.DOF(FEA.x == 0 & FEA.y == 0,:)'; % corner node is clamped
        
        foo = FEA.DOF(FEA.x <= wing.root_chord/4 & FEA.x > 0 & FEA.y == 0,:);
        FEA.bc = [FEA.bc;foo(:)]; % portion of root is clamped
        
        foo = FEA.DOF(abs(FEA.x-ppval(wing.LE_pp,FEA.y)) < 1E-4 & FEA.y > 0 & FEA.y <= wing.root_chord/2,:);
        FEA.bc = [FEA.bc;foo(:)]; % portion of LE is clamped
        
        foo = FEA.DOF(setxor(1:FEA.N_nodes,FEA.beam_elements(:)),4:6)';
        FEA.bc = [FEA.bc;foo(:)]; %membranes have no rotational DOF
        
        FEA.bc = unique(FEA.bc);
        FEA.fdof = 1:FEA.GlobalDOF; FEA.fdof(FEA.bc(:,1)) = [];
        clear foo
        
        %% Stiffness matrix and conversion of pressures to forces:
        
        [FEA.K,interp.F_dP] = stiffness_matrix(FEA.x,FEA.y,FEA.z,FEA.trielements,FEA.beam_elements,FEA.membrane_asmb,FEA.beam_asmb,FEA.beam_ID,FEA.GlobalDOF,FEA.fdof,mat);
        
        %% Mass matrix:
        
        FEA.M = mass_matrix(FEA.x,FEA.y,FEA.z,FEA.trielements,FEA.beam_elements,FEA.membrane_asmb,FEA.beam_asmb,FEA.beam_ID,FEA.GlobalDOF,FEA.fdof,mat);
        
        %% Natural modes:
        
        options.disp = 0;
        [interp.modes,natural] = eigs(FEA.K,FEA.M,FEA.N_modes,'sm',options);
        natural = sqrt(diag(natural));
        [natural,ID] = sort(natural);
        
        interp.modes = interp.modes(:,ID);
        foo = interp.modes'*FEA.M*interp.modes;
        for j = 1:FEA.N_modes
            interp.modes(:,j) = interp.modes(:,j)/sqrt(foo(j,j));
        end
        FEA.naturals = natural;
        clear ID foo i j natural
        
        FEA.Mr = eye(FEA.N_modes);
        FEA.Cr = FEA.Mr*mat.damping;
        FEA.Kr = diag(FEA.naturals.^2);
        interp.F_dP = interp.modes'*interp.F_dP;
        
        %% Inertial forces:
        
        if kin.periodic
            [FEA.F_iner,loads.CP_inertial] = inertial_loads(FEA.x+wing.x_offset,FEA.y+wing.y_offset,FEA.z,FEA.trielements,FEA.beam_elements,FEA.membrane_asmb,FEA.beam_asmb,FEA.beam_ID,FEA.GlobalDOF,FEA.fdof,mat,...
                kin.zpp(1:time.N_per)',kin.f(1:time.N_per)',kin.fp(1:time.N_per)',kin.fpp(1:time.N_per)',kin.th(1:time.N_per)',kin.thp(1:time.N_per)',kin.thpp(1:time.N_per)');
            FEA.F_iner = interp.modes'*FEA.F_iner;
            FEA.F_iner = [repmat(FEA.F_iner,1,time.N_cycles),FEA.F_iner(:,1)];
            loads.CP_inertial = [repmat(loads.CP_inertial,1,time.N_cycles),loads.CP_inertial(1)];
        else
            [FEA.F_iner,loads.CP_inertial] = inertial_loads(FEA.x+wing.x_offset,FEA.y+wing.y_offset,FEA.z,FEA.trielements,FEA.beam_elements,FEA.membrane_asmb,FEA.beam_asmb,FEA.beam_ID,FEA.GlobalDOF,FEA.fdof,mat,...
                kin.zpp',kin.f',kin.fp',kin.fpp',kin.th',kin.thp',kin.thpp');
            FEA.F_iner = interp.modes'*FEA.F_iner;
        end
        
        loads.CP_inertial = loads.CP_inertial'./(sum(2*gauss.b*wing.length/gauss.N_stations)*.5*setup.rho*setup.U.^3);
        
        %% Convert pressures at Gauss points into uniform pressures over elements:
        
        interp.dP_convert = pressure_interp(gauss.x,gauss.y,gauss.b,gauss.LE,gauss.N_gauss,gauss.N_stations,wing.LE_knots,wing.TE_knots,FEA.x,FEA.y,FEA.trielements);
        
        %% Convert FEA solution into displacements at gauss points:
        
        interp.u_convert = disp_interp(gauss.x,gauss.y,gauss.b,gauss.LE,gauss.N_gauss,gauss.N_stations,FEA.x,FEA.y,FEA.trielements,FEA.bc,FEA.fdof);
        interp.u_convert = interp.u_convert*interp.modes;
        
        %% Initial conditions:
        
        solve.X = zeros(inflow.N*gauss.N_stations+2*FEA.N_modes,1);
        
        %% Allocate:
        
        FEA.mode_history = zeros(length(time.t),5); % store the first 5 modal amplitudes
        loads.CL = zeros(length(time.t),1); % lift coefficient
        loads.CD = zeros(length(time.t),1); % drag coefficient
        loads.CP_aero = zeros(length(time.t),1); % aerodynamic power
        loads.CP_KER = zeros(length(time.t),1); % kinetic energy rate
        loads.CP_SER = zeros(length(time.t),1); % strain energy rate
        loads.CP_total = zeros(length(time.t),1); % total power coefficient (inertial + aero + KER + SER)
        
        %% Pre-compute some matrices:
        
        interp.T_u = interp.T_hn_h*interp.u_convert;
        interp.T_zo = interp.T_hn_h*gauss.zo;
        interp.K_T_zo = aero.K*interp.T_hn_h*gauss.zo;
        interp.W_K_T_zo = aero.W*interp.K_T_zo;
        interp.K_T_u = aero.K*interp.T_hn_h*interp.u_convert;
        interp.W_K_T_u = aero.W*interp.K_T_u;
        interp.F_map = 2*setup.rho*interp.F_dP*interp.dP_convert*loads.correct_total*interp.P_map;
        interp.P_map2 = 2*setup.rho*loads.correct_total*interp.P_map;
        
        inflow.Rupp = inflow.Rw*interp.T_hn_h*interp.u_convert;
        FEA.M_eff = FEA.Mr - 2*setup.rho*interp.F_dP*interp.dP_convert*loads.correct_total*interp.P_map*aero.W*interp.T_hn_h*interp.u_convert;
        
        %% Time marching:
        
        if plot_decide
            figure
            set(gcf,'position',[98 160 1324 938])
        end
        
        for i = 1:length(time.t)
            
            %% Transformation matrix:
            
            kin.T = [cos(kin.th(i)),0,sin(kin.th(i));0,1,0;-sin(kin.th(i)),0,cos(kin.th(i))]*[1,0,0;0,cos(kin.f(i)),-sin(kin.f(i));0,sin(kin.f(i)),cos(kin.f(i))];
            
            %% Compute aeroelastic matrices:
            
            aero.uo_mult = kron(diag(aero.uo(i,:)),speye(aero.N_h+1));
            aero.uo_dot_mult = kron(diag(aero.uo_dot(i,:)),speye(aero.N_h+1));
            
            inflow.B = kron(diag(aero.uo(i,:)./gauss.b),speye(inflow.N));
            inflow.R = inflow.Rw*(aero.vn_dot(:,i)+aero.uo_dot_mult*interp.K_T_zo);
            inflow.Ru = inflow.Rw*aero.uo_dot_mult*interp.K_T_u;
            inflow.Rup = inflow.Rw*aero.uo_mult*interp.K_T_u;
            FEA.K_eff = FEA.Kr - interp.F_map*(aero.uo_mult*aero.uo_mult*interp.K_T_u + aero.uo_dot_mult*interp.W_K_T_u);
            FEA.C_eff = FEA.Cr - interp.F_map*aero.uo_mult*(interp.T_u + interp.W_K_T_u);
            FEA.F_aero = interp.F_map*(aero.uo_mult*aero.vn(:,i) + aero.W*aero.vn_dot(:,i) + aero.uo_mult*aero.uo_mult*interp.K_T_zo + aero.uo_dot_mult*interp.W_K_T_zo);
            FEA.F_lambda = -interp.F_map*aero.uo_mult*inflow.Q;
            
            %% Assemble:
            
            solve.A = [inflow.A,zeros(inflow.N*gauss.N_stations,FEA.N_modes),-inflow.Rupp;zeros(FEA.N_modes,inflow.N*gauss.N_stations),eye(FEA.N_modes),zeros(FEA.N_modes);zeros(FEA.N_modes,inflow.N*gauss.N_stations),zeros(FEA.N_modes),FEA.M_eff];
            solve.B = [inflow.B,-inflow.Ru,-inflow.Rup;zeros(FEA.N_modes,inflow.N*gauss.N_stations),zeros(FEA.N_modes),-eye(FEA.N_modes);-FEA.F_lambda,FEA.K_eff,FEA.C_eff];
            solve.R = [inflow.R;zeros(FEA.N_modes,1);FEA.F_iner(:,i)+FEA.F_aero];
            
            %% Integrate:
            
            if i == 1
                solve.Xp = solve.A\(solve.R-solve.B*solve.X);
                solve.G = ((1/solve.beta/time.t_step)*solve.A + solve.alpha*solve.B)*solve.X + (1/solve.beta-1)*solve.A*solve.Xp - solve.alpha*solve.R;
            else
                solve.X_old = solve.X;
                solve.Xp_old = solve.Xp;
                solve.X = ((1/solve.beta/time.t_step)*solve.A + (1+solve.alpha)*solve.B)\((1+solve.alpha)*solve.R + solve.G);
                solve.Xp = (1/solve.beta/time.t_step)*(solve.X-solve.X_old)+(1-1/solve.beta)*solve.Xp_old;
                solve.G = ((1/solve.beta/time.t_step)*solve.A + solve.alpha*solve.B)*solve.X + (1/solve.beta-1)*solve.A*solve.Xp - solve.alpha*solve.R;
            end
            
            %% Inflow terms:
            
            inflow.lambda = reshape(solve.X(1:inflow.N*gauss.N_stations),inflow.N,gauss.N_stations);
            aero.lambda_o = .5*inflow.b'*inflow.lambda;
            aero.lambda = zeros((aero.N_h+1)*gauss.N_stations,1); aero.lambda(1:(aero.N_h+1):size(aero.vn,1),:) = aero.lambda_o';
            
            %% Structural terms:
            
            FEA.eta = solve.X(inflow.N*gauss.N_stations+1:inflow.N*gauss.N_stations+FEA.N_modes);
            FEA.eta_p = solve.X(inflow.N*gauss.N_stations+FEA.N_modes+1:end);
            FEA.eta_pp = solve.Xp(inflow.N*gauss.N_stations+FEA.N_modes+1:end);
            FEA.answer = zeros(FEA.GlobalDOF,1); FEA.answer(FEA.fdof) = interp.modes*FEA.eta; FEA.answer = reshape(FEA.answer,6,FEA.N_nodes)';
            FEA.mode_history(i,:) = FEA.eta(1:5)';
            
            %% Airloads:
            
            aero.hn = interp.T_zo + interp.T_u*FEA.eta;
            aero.hn_dot = interp.T_u*FEA.eta_p;
            aero.hn_dot_dot = interp.T_u*FEA.eta_pp;
            aero.wn = aero.vn(:,i) + aero.hn_dot + aero.uo_mult*aero.K*aero.hn;
            aero.wn2 = aero.vn(:,i) + aero.hn_dot; aero.wn2(1:(aero.N_h+1):size(aero.vn,1)) = aero.wn2(1:(aero.N_h+1):size(aero.vn,1)) - .5*setup.U(i)*(sin(kin.f(i)+kin.th(i))-sin(kin.f(i)-kin.th(i)));
            
            aero.taun = aero.uo_mult*(aero.vn(:,i)-aero.lambda) + aero.W*aero.vn_dot(:,i) + ...
                (aero.uo_mult*aero.uo_mult)*(interp.K_T_zo + interp.K_T_u*FEA.eta) + aero.uo_dot_mult*(interp.W_K_T_zo + interp.W_K_T_u*FEA.eta) + ...
                aero.uo_mult*(interp.T_u + interp.W_K_T_u)*FEA.eta_p + aero.W*aero.hn_dot_dot;
            aero.dP = interp.P_map2*aero.taun;
            gauss.hx = reshape(interp.hx_map*aero.hn,gauss.N_gauss,gauss.N_stations);
            gauss.dP = reshape(aero.dP,gauss.N_gauss,gauss.N_stations);
            
            %% Final forces:
            
            loads.alpha = reshape(-interp.h_map*aero.hn,gauss.N_gauss,gauss.N_stations);
            loads.alpha = atan2(loads.alpha(end,:)-loads.alpha(1,:),2*gauss.b) + atan2(aero.vo(i,:),aero.uo(i,:));
            loads.fz = gauss.b.*sum(repmat(gauss.weight,1,gauss.N_stations).*gauss.dP.*repmat(sin(gauss.theta),1,gauss.N_stations)) + ...
                setup.rho*setup.U(i)*setup.U(i)*loads.correct_station.*gauss.b.*(loads.cd_o*cos(loads.alpha).^2+loads.cd_pi2*sin(loads.alpha).^2).*aero.vo(i,:)./(sqrt(aero.vo(i,:).^2+aero.uo(i,:).^2));
            loads.fx = gauss.b.*sum(repmat(gauss.weight,1,gauss.N_stations).*gauss.dP.*gauss.hx.*repmat(sin(gauss.theta),1,gauss.N_stations)) + ...
                -2*pi*setup.rho*loads.correct_station.*gauss.b.*(aero.wn(1:aero.N_h+1:end)'-aero.lambda_o).^2 + ...
                setup.rho*setup.U(i)*setup.U(i)*loads.correct_station.*gauss.b.*(loads.cd_o*cos(loads.alpha).^2+loads.cd_pi2*sin(loads.alpha).^2).*aero.uo(i,:)./(sqrt(aero.vo(i,:).^2+aero.uo(i,:).^2));
            loads.temp = kin.T*[loads.fx;zeros(1,gauss.N_stations);loads.fz];
            loads.Fx = loads.temp(1,:); loads.Fy = loads.temp(2,:); loads.Fz = loads.temp(3,:);
            
            loads.cl = loads.Fz./gauss.b/setup.rho/setup.U(i)/setup.U(i);
            loads.cd = loads.Fx./gauss.b/setup.rho/setup.U(i)/setup.U(i);
            loads.CL(i,1) = (wing.length/gauss.N_stations)*sum(loads.cl)/wing.length;
            loads.CD(i,1) = (wing.length/gauss.N_stations)*sum(loads.cd)/wing.length;
            
            %% Power:
            
            loads.cp_aero = sum(repmat(gauss.weight,1,gauss.N_stations).*reshape(interp.h_map*aero.wn2,gauss.N_gauss,gauss.N_stations).*gauss.dP.*repmat(sin(gauss.theta),1,gauss.N_stations))/setup.rho/setup.U(i)/setup.U(i)/setup.U(i);
            loads.CP_aero(i,1) = (wing.length/gauss.N_stations)*sum(loads.cp_aero)/wing.length;
            loads.CP_KER(i,1) = (FEA.eta_p'*(FEA.Mr*FEA.eta_pp + FEA.Cr*FEA.eta_p))/(sum(2*gauss.b*wing.length/gauss.N_stations)*.5*setup.rho*setup.U(i)^3);
            loads.CP_SER(i,1) = (FEA.eta_p'*FEA.Kr*FEA.eta)/(sum(2*gauss.b*wing.length/gauss.N_stations)*.5*setup.rho*setup.U(i)^3);
            loads.CP_total(i,1) = loads.CP_inertial(i,1) + loads.CP_KER(i,1) + loads.CP_SER(i,1) + loads.CP_aero(i,1);
            
            %% Plotting:
            
            if plot_decide
                if mod(i,time.skip) == 0 || i == 1 || i == length(time.t)
                    
                    subplot(2,3,1),plot3(wing.LE_knots,wing.y_stations,wing.y_stations*0,'b.',wing.TE_knots,wing.y_stations,wing.y_stations*0,'b.','markersize',15)
                    hold on
                    plot3(wing.x(1,:),wing.y(1,:),0*wing.x(1,:),'b-',wing.x(end,:),wing.y(end,:),0*wing.x(end,:),'b-','linewidth',2)
                    surf(wing.x,wing.y,wing.z,'facecolor',[.5 .5 .5],'edgecolor','none','facealpha',.5)
                    plot3(gauss.x+reshape(repmat(gauss.b,gauss.N_gauss,1),1,gauss.N_gauss*gauss.N_stations)'+reshape(repmat(gauss.LE,gauss.N_gauss,1),1,gauss.N_gauss*gauss.N_stations)',reshape(repmat(gauss.y,gauss.N_gauss,1),1,gauss.N_gauss*gauss.N_stations)',-gauss.zo,'r.')
                    axis equal, axis tight
                    axis([min(min(wing.x)) max(max(wing.x)) 0 max(max(wing.y)) min(min(wing.z))-.01 max(max(wing.z))+.01])
                    xlabel('x')
                    ylabel('y')
                    view([30 34])
                    
                    
                    for j=1:connectivity.edges.Nedges
                        if connectivity.edges.type(j) == 'I'
                            subplot(4,3,2),plot([connectivity.vertices.coords(connectivity.edges.vertices(j,1),1),connectivity.vertices.coords(connectivity.edges.vertices(j,2),1)],[connectivity.vertices.coords(connectivity.edges.vertices(j,1),2),connectivity.vertices.coords(connectivity.edges.vertices(j,2),2)],'r-','MarkerEdgeColor','k','LineWidth',.5); hold on;
                        else
                            subplot(4,3,2),plot([connectivity.vertices.coords(connectivity.edges.vertices(j,1),1),connectivity.vertices.coords(connectivity.edges.vertices(j,2),1)],[connectivity.vertices.coords(connectivity.edges.vertices(j,1),2),connectivity.vertices.coords(connectivity.edges.vertices(j,2),2)],'k-','MarkerEdgeColor','k','LineWidth',2); hold on;
                        end
                    end
                    axis equal; axis([-1.01 1.01 -1.01 1.01]); axis off
                    
                    
                    subplot(4,3,5),trisurf(FEA.trielements,FEA.x,FEA.y,FEA.z,'facecolor','none','edgecolor','k')
                    hold on
                    for j = 1:length(FEA.beam_elements)
                        if FEA.beam_ID(j,:) == 'W' || FEA.beam_ID(j,:) == 'N' || FEA.beam_ID(j,:) == 'E'
                            plot3([FEA.x(FEA.beam_elements(j,1)),FEA.x(FEA.beam_elements(j,2))],[FEA.y(FEA.beam_elements(j,1)),FEA.y(FEA.beam_elements(j,2))],[FEA.z(FEA.beam_elements(j,1)),FEA.z(FEA.beam_elements(j,2))],'k-','linewidth',2)
                        elseif FEA.beam_ID(j,:) == 'S'
                            plot3([FEA.x(FEA.beam_elements(j,1)),FEA.x(FEA.beam_elements(j,2))],[FEA.y(FEA.beam_elements(j,1)),FEA.y(FEA.beam_elements(j,2))],[FEA.z(FEA.beam_elements(j,1)),FEA.z(FEA.beam_elements(j,2))],'b-','linewidth',2)
                        else
                            plot3([FEA.x(FEA.beam_elements(j,1)),FEA.x(FEA.beam_elements(j,2))],[FEA.y(FEA.beam_elements(j,1)),FEA.y(FEA.beam_elements(j,2))],[FEA.z(FEA.beam_elements(j,1)),FEA.z(FEA.beam_elements(j,2))],'r-','linewidth',2)
                        end
                    end
                    axis equal, axis tight
                    axis([min(min(wing.x)) max(max(wing.x)) 0 max(max(wing.y)) min(min(wing.z))-.01 max(max(wing.z))+.01])
                    view([90 90])
                    grid off
                    
                    
                    subplot(4,3,3),surf(wing.x,wing.y,wing.z,'edgecolor','none','facecolor',[.5 .5 .5],'facealpha',.5)
                    foo = reshape(gauss.x+reshape(repmat(gauss.b,gauss.N_gauss,1),1,gauss.N_gauss*gauss.N_stations)'+reshape(repmat(gauss.LE,gauss.N_gauss,1),1,gauss.N_gauss*gauss.N_stations)',gauss.N_gauss,gauss.N_stations);
                    hold on
                    for j = 1:gauss.N_stations
                        plot3(foo(1:end-2,j),repmat(gauss.y(j),gauss.N_gauss-2,1),gauss.dP(1:end-2,j)/(.5*setup.rho*setup.U(i)*setup.U(i))/200,'k-','linewidth',1.5)
                    end
                    clear temp
                    axis equal, axis tight
                    axis([0 wing.root_chord 0 wing.length -.1 .1])
                    zlabel('\DeltaC_p/200')
                    view([30 34])
                    set(gca,'position',[0.6677    0.7452    0.2553    0.2420])
                    
                    
                    subplot(4,3,6),plot(gauss.y/wing.length,loads.cl,'k.-',gauss.y/wing.length,loads.cd,'r.-')
                    xlabel y/L
                    ylabel('c_l, c_d')
                    axis([0 1 -3 3])
                    
                    
                    foo = kin.T*[FEA.x'+wing.x_offset;FEA.y'+wing.y_offset;FEA.z'];
                    wing.Xo = foo(1,:)'; wing.Yo = foo(2,:)'; wing.Zo = foo(3,:)' + kin.z(i);
                    foo = kin.T*[FEA.x'+FEA.answer(:,1)'+wing.x_offset;FEA.y'+FEA.answer(:,2)'+wing.y_offset;FEA.z'+FEA.answer(:,3)'];
                    wing.X = foo(1,:)'; wing.Y = foo(2,:)'; wing.Z = foo(3,:)' + kin.z(i);
                    subplot(2,2,3),trisurf(FEA.trielements,wing.Xo,wing.Yo,wing.Zo,FEA.answer(:,3),'edgecolor','none','facecolor',[.5 .5 .5],'facealpha',0.5)
                    hold on
                    trisurf(FEA.trielements,wing.X,wing.Y,wing.Z,FEA.answer(:,3),'edgecolor','none','facecolor','interp')
                    axis equal
                    axis([-.1 .2 0 .5 -.3 .3])
                    xlabel('x')
                    ylabel('y')
                    view([30 34])
                    clb = colorbar('vert');
                    set(gca,'position',[0.0181    0.0330    0.3927    0.5650])
                    set(clb,'position',[0.0219    0.1321    0.0102    0.2772])
                    
                    
                    subplot(10,10,60),plot(time.t(1:i)*kin.omega/2/pi,loads.CL(1:i),'k-')
                    xlabel t/T
                    ylabel C_L
                    set(gca,'position',[0.4108    0.3291    0.2134    0.1577])
                    
                    
                    subplot(4,3,9),plot(time.t(1:i)*kin.omega/2/pi,loads.CP_aero(1:i),'r-',...
                        time.t(1:i)*kin.omega/2/pi,loads.CP_KER(1:i),'b-',...
                        time.t(1:i)*kin.omega/2/pi,loads.CP_inertial(1:i),'g-',...
                        time.t(1:i)*kin.omega/2/pi,loads.CP_SER(1:i),'m-')
                    hold on
                    plot(time.t(1:i)*kin.omega/2/pi,loads.CP_total(1:i),'k-','linewidth',2)
                    xlabel t/T
                    ylabel C_P
                    
                    
                    subplot(10,10,90),plot(time.t(1:i)*kin.omega/2/pi,loads.CD(1:i),'k-')
                    xlabel t/T
                    ylabel C_D
                    set(gca,'position',[0.4108    0.1100    0.2134    0.1577])
                    
                    
                    subplot(4,3,12),plot(time.t(1:i)*kin.omega/2/pi,FEA.mode_history(1:i,:))
                    xlabel t/T
                    ylabel \eta
                    
                    
                    pause(0.1)
                    if i < length(time.t)
                        clf
                    end
                    
                end
            end
            
        end
        
        %% Efficiency:
        
        loads.eta_total = -mean(loads.CD(length(time.t)-time.N_per:length(time.t)))/mean(loads.CP_total(length(time.t)-time.N_per:length(time.t)));
        loads.eta_aero = -mean(loads.CD(length(time.t)-time.N_per:length(time.t)))/mean(loads.CP_aero(length(time.t)-time.N_per:length(time.t)));
        
        %% Final objectives (cycle-averaged drag, power, and lift)
        
        g = [mean(loads.CD(length(time.t)-time.N_per:length(time.t))),mean(loads.CP_total(length(time.t)-time.N_per:length(time.t))),-mean(loads.CL(length(time.t)-time.N_per:length(time.t)))];
        
    else
        
        g = g_fail;
        
    end
    
else
    
    g = g_fail;
    
end

