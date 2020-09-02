
% Title: Enforcing Safety at Runtime for Systems with Disturbances
% Submitted to 2020 IEEE 59th Conference on Decision and Control (CDC)
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 3/28/2020
% Description:  This script generates the figures for the case study.
%               An ASIF is created to assure a nondeterministic system.

clf; clear all; clc;


% -------------------------------------------------------------------------
% User Inputs: Parameters for Dynamics
% -------------------------------------------------------------------------

% Description of system:
%    x1 = velocity of cart 1
%    x2 = velocity of cart 2
%    x3 = velocity of cart 3
%    x4 = position cart 1 - position cart 2
%    x5 = position cart 2 - position cart 3

global W beta
n = 5; % 5 states (3 velocities 2 displacements)
m = 2; % 2 inputs (2 nonlinear springs)
       % 3 disturbances (effect the velcoity states)
     
W = [-.1, .1; ...
     -.1, .1; ...
     -.1, .1];      % Disturbance Bounds

beta = -1;             % Friction term (always negitive)


% -------------------------------------------------------------------------
% User Inputs: Parameters for Backup Controller
% -------------------------------------------------------------------------
global Springs

Springs = [2, .5];  % contains spring parameters for backup controller
                    % [Hooke's Const, Saturation] (always positive)
                      
% -------------------------------------------------------------------------
% Construct Graph Structure
% -------------------------------------------------------------------------
global N K D
N = n - m; % 3 carts
K = m;     % 2 edges for comunication
D = [-1,  0; ...
      1, -1; ...
      0,  1]; % connectivity matrix
clear n m

% -------------------------------------------------------------------------
% User Inputs: Video Parameters
% -------------------------------------------------------------------------
global create_video myVideo
create_video = 0;

% -------------------------------------------------------------------------
% User Inputs: Simulation Parameters
% -------------------------------------------------------------------------
global Tb x0 T dt
x0 = [1/2; 3/7 ; 5/4; 1/4; 1/2];  % Initial condition
Tb = 1;                           % Backup time
T  = 4;                           % Simulation Time
dt = .01;                         % Timestep for Simulation
T_sim = 0:dt:T;                   % Vector of Simulation Times


% -------------------------------------------------------------------------
% Get Random Continuous Disturbance Signal from Gaussian Process
% -------------------------------------------------------------------------
% Generate random set of velocity inputs to the unicylcye
n=length(T_sim);
param=1;
K_ss=kernel(T_sim,T_sim, param);
CHOL =chol(K_ss+1e-12*eye(n),'lower');
% Generate 3 random deisturbance signals from GP
for i = 1:3
    gp_signal = (CHOL*normrnd(0, 1, length(T_sim),1))';
    const = (W(i, 2)- W(i, 1))/(max(gp_signal) - min(gp_signal));
    w_signal(i, :) = const*(gp_signal - min(gp_signal))  + W(i, 1);
end

% -------------------------------------------------------------------------
% Get Safe Set and Verify that Sb+(Tb) does not Intersect Unsafe Set 
% -------------------------------------------------------------------------

global P L
[P, L] = get_Safe_set();            % Find Parameters for Elipitcal Safe Set 

% Overapproximate Safe Set with rectangle
[S]    = over_approx_Sb();          % S is a rectiable overapproximating Sb
S = [ floor(S(:, 1)) + floor((S(:, 1) - floor(S(:, 1)))/0.01) * 0.01, ... % round down a bit
      floor(S(:, 2)) + ceil( (S(:, 2) - floor(S(:, 2)))/0.01) * 0.01];    % round up a bit


[~, QQ] = PhiE2(Tb, S(:, 1), S(:, 2));
% Approx is a hyper rectangular  overapproximation of the Tb sectond
% backward reachable set of S under the backup dynamics.  Thus if Approx
% does not intersect the unsafe set, then the Tb second basin of attracton
% of Sb under the backup dynamics is safe under the backup controller
Approx = [QQ(1:5, end), QQ(6:10, end)];
fprintf('Tb = 1 verified to be safe backup time\n');

% -------------------------------------------------------------------------
% Setup plots
% -------------------------------------------------------------------------

% this plot shows the vehicle trajectories
figure(1); clf; 

ploot(1)
hold on; 
axis([-.25, 2.25, -.1, .25])
set(gca,'XTick',0:.5:2, 'YTick', [])
plot([-.5, 2.5], [-.025, -.025], 'k', 'LineWidth', 2)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
Leg = legend();
set(Leg,'visible','off');
xlabel({'Distance', '', ''},'Interpreter','latex')

ploot(2)
hold on; grid on
xlabel({'$x_4$', ''},'Interpreter','latex')
ylabel('$x_5$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
Leg = legend();
set(Leg,'visible','off');
axis([0, 2, -.25, .75])
xticks(0:.5:2)
yticks(-.25:.25:.75)

% this plot shows the applied input vs time
% Input 1
ploot(3); hold on; grid on
xlabel('$t$','Interpreter','latex')
ylabel('$u_{\rm d},\, u^{\rm ASIF}$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
Leg = legend();
set(Leg,'visible','off');
axis([0, T, -.35, .35])
yticks(-.3:.3:.3)

% Input 2
ploot(4); hold on; grid on
xlabel('$t$','Interpreter','latex')
ylabel('$u_{\rm d},\, u^{\rm ASIF}$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
Leg = legend();
set(Leg,'visible','off');
axis([0, T, -.3, .3])
yticks(-.3:.3:.3)


% -------------------------------------------------------------------------
% Simulate carts moving with controller
% -------------------------------------------------------------------------

% Initialise
STATE_TRAJ = x0;
STATE_TRAJ_nom = STATE_TRAJ;
U_Sim = [];            

% Set up meshgrid for plottng elliptical safe set
np = 150;
[XX1,XX2] = meshgrid(linspace(-2,2,np),linspace(-2,2,np));


% the stopwatch function calculates the time for each optimazation
global comp_times
comp_times = [];    % this variable stores the time of each optimazation


% set up cart picture
CART_SCALE = .1;% For plotting graphic
CART = MakeCart(CART_SCALE);

plot_holder  = []; % stores plot variables
plot_holder2 = []; % stores plot variables
plot_holder3 = []; % stores plot variables

fprintf('Begin Simulation: \n');
sz = T/dt;          % number of time steps for simulation
video_maker(1);     % set up video

cart1_pos = 0;
for t = 1:sz
    fprintf(['Timestep: ', int2str(t),'/ ', int2str(sz), ' \n']); % print iteration
    mean(comp_times)
    %stopwatch()
    
    w = w_signal(:, t);  % choose disturbances
    
    u_des = [-.3*sin((1/2)*t*(2*pi/sz)); ...
             .2*cos(t*(2*pi/sz))]; % desired input
    
    [u_app, xmax, Phi] = QP(STATE_TRAJ(:, t), u_des);    % run asif QP
    
    % store applied input and next state
    STATE_TRAJ(:, t+1) = Simulate_Dynamics(STATE_TRAJ(:, t), u_app, w);
    STATE_TRAJ_nom(:, t+1) = Simulate_Dynamics(STATE_TRAJ_nom(:, t), u_des, w);
    
    U_Sim(:, t) = [u_des; u_app];
    
    %%%%%%%%%%
    % Plot
    %%%%%%%%%%%
    if t~=1
        delete(plot_holder)
    end
    plot_holder = [];
    
    ploot(1)
    cart1_pos = cart1_pos + dt*STATE_TRAJ(1, t);
    cart2_pos = cart1_pos + STATE_TRAJ(4, t);
    cart3_pos = cart2_pos + STATE_TRAJ(5, t);
    
    cart1 = plot(cart1_pos+CART(1, :), CART(2, :), 'r', ...
        'LineWidth', 2);
    cart2 = plot(cart2_pos+CART(1, :), CART(2, :), 'b', ...
        'Linewidth', 2);
    cart3 = plot(cart3_pos+CART(1, :), CART(2, :), 'Color', [0, .8, .1], ...
        'Linewidth', 2);
    
    plot_holder = [cart1, cart2, cart3];
    
    ploot(2)
    % For the current velocities, project the Safe Set to the x4-x5 plane 
    corners   = get_corners(xmax(1:5),xmax(6:10));
    hs_corners=cellfun(@hsam,num2cell(corners,1));
    [~,worst_corner]=min(hs_corners); %best case safety
    xrect     = corners(:,worst_corner);
    tempvect  = xrect;
    for c1=1:np
        for c2=1:np
            tempvect(4:5)=[XX1(c1,c2),XX2(c1,c2)];
            ZZ(c1,c2)=hsam(tempvect);
        end
    end
    clear c1 c2
    
    % Plot Safe Set
    [cc, hh] = contour(XX1, XX2, ZZ, ...
                    [0 0], 'r', ...
                    'HandleVisibility', 'off'); % Overlay contour line
    plot_holder(end+1) = patch(cc(1,2:end), cc(2, 2:end), [0, 1, 0], ...
                    'EdgeColor', [0, 0, 0], ...
                    'FaceAlpha', .1, ...
                    'HandleVisibility', 'off');
    delete(hh);
                      
    % Best Case Approximation of Reachable Set (Rectangle)           
    plot_holder(end+1:end+2) = plot_rect(xmax(4:5), xmax(9:10));

    % Plot total reachable set approximation on (0 - Tb)
    plot_holder(end+1) = plot([flip(Phi(4, :)), Phi( 9, :)], ...
                              [flip(Phi(5, :)), Phi(10, :)], ...
                              'Color', [1, 0, 0, .4], ...
                              'LineWidth', 2);
                          
    % Plot corners of reachble set and worst case corner                      
    plot_holder(end+1) = scatter(xmax([4, 9]), xmax([5, 10]), 'k', 'filled');
    plot_holder(end+1) = scatter(xrect(4), xrect(5), ...
                            'MarkerFaceColor', 'r', ...
                            'MarkerEdgeColor', 'k');
                        
    % Nominal Trajectory
    plot_holder(end+1) = plot(STATE_TRAJ_nom(4, :), STATE_TRAJ_nom(5, :), ...
                            'm', 'LineWidth', 2);
    plot_holder(end+1) = scatter(STATE_TRAJ_nom(4, end), STATE_TRAJ_nom(5, end), ...
                            'm','filled');
                        
    % ASIF Trajectory
    plot_holder(end+1) = plot(STATE_TRAJ(4, :), STATE_TRAJ(5, :), ...
                            'b', 'LineWidth', 2);
    plot_holder(end+1) = scatter([STATE_TRAJ(4, 1), STATE_TRAJ(4, end)], ...
                                 [STATE_TRAJ(5, 1), STATE_TRAJ(5, end)], ...
                            'b', 'filled');
                        
    drawnow;
    %%%%%%%
    
    ploot(3)
    if t~=1
        delete(plot_holder2);
        plot_holder2 = [];
    end
    plot_holder2(1) = plot(T_sim(1, 1:t), U_Sim(1, :), 'r', 'LineWidth', 2);
    plot_holder2(2) = plot(T_sim(1, 1:t), U_Sim(3, :), 'b', 'LineWidth', 2);
    drawnow;
    
    ploot(4)
    if t~=1
        delete(plot_holder3);
        plot_holder3 = [];
    end
    plot_holder3(1) = plot(T_sim(1, 1:t), U_Sim(2, :), 'r', 'LineWidth', 2);
    plot_holder3(2) = plot(T_sim(1, 1:t), U_Sim(4, :), 'b', 'LineWidth', 2);
    drawnow;
    
    video_maker(2); % Record Frame
end
video_maker(3); % Close Video

fprintf('Simulation completed.  Average solver time was: \n');
mean(comp_times)

%------------------------------------------------------------
% 3: Functions
%------------------------------------------------------------

% -------------
%  control
% -------------
function [u_app, xmax, Phi] = QP(x, u_des)
    global W beta Tb N K D 
    global comp_times % for storing cpu times for optimazation
    
    % gets every possible combination of worst case disturbances
    w = get_corners(W(:, 1), W(:, 2));
    
    % evaluate the state flow in the embedding space
    % evaluate the jacobian of the state flow
    [Tout, Phi] = PhiE(Tb, x, x);
    
    hs_soft = cellfun(@hsamsoftmin,num2cell(Phi,1)); %when this is positive, all corners are within the safe backup region
    
    % maxtimeind = index of max time
    [Psi, maxtimeind] = max(hs_soft); %best case safety

    DPsi = hsamsoftmingrad(Phi(:,maxtimeind))*DPhiE(Tout(maxtimeind),x);

    %%%%%
    xmax      = Phi(:,maxtimeind);
    corners   = get_corners(xmax(1:5),xmax(6:10));
    hs_corners=cellfun(@hsam,num2cell(corners,1));
    [~,worst_corner]=min(hs_corners); %best case safety
    xrect     = corners(:,worst_corner);
    np        = 150;
    [XX1,XX2] = meshgrid(linspace(-2,2,np),linspace(-2,2,np));
    tempvect  = xrect;
    for c1=1:np
        for c2=1:np
            tempvect(4:5)=[XX1(c1,c2),XX2(c1,c2)];
            ZZ(c1,c2)=hsam(tempvect);
        end
    end
    
    %%%%
    
    A = [beta*eye(N), zeros(N, K); D', zeros(K)];
    B = [D; zeros(K)];
    C = [eye(N); zeros(K, N)];

    cvx_begin quiet
        %Optimize over u
        variable u(K, 1)
        
        % QP Constraints
        for i = 1:2^N
            dxdt = A*x - B*u + C*w(:, i);
            DPsi*dxdt >=- 1000*Psi^3; % constraint on embedding x
        end
        minimize( norm(u - u_des, 2) );
    cvx_end;
    
    comp_times = [comp_times, cvx_cputime]; % save computation time
    
    if isequal(cvx_status, 'Failed') || isequal(cvx_status, 'Infeasible')
        error('QP Error')
    end
    
    u_app = u;
end

function u = ub(x)
    global K N Springs

    holder =zeros(K, 1);
    for i = 1:K
        holder(i, 1) = Springs(1)*tanh(Springs(2)*x(N+i));
    end
    u = holder;
end

function out = Simulate_Dynamics(xt, u, w)
    global N K D beta dt
    A = [beta*eye(N), zeros(N, K); D', zeros(K)];
    B = [D; zeros(K)];
    C = [eye(N); zeros(K, N)];
    
    f = 0;
    for tau = linspace(0, dt, 201)
        f = f + (dt/200)*expm(A*(dt - tau));
    end
    % the control input is constant over dt so this makes sense
    out = expm(A*dt)*xt + f*(-B*u + C*w); % next state
end

% -------------
%  estimation
% -------------

function out = decomp(x, xh, w)
    % Decompisition function for backup dynamics
    global beta N K D 
    
    Dplus = (D > 0).*D;  % Get positive part of D
    Dminus = (D < 0).*D; % Get negitive part of D

    A1 = [beta.*eye(N), zeros(N, K); Dplus', zeros(K)];
    A2 = [zeros(N), zeros(N, K); Dminus', zeros(K)];
    B1 = [-Dminus; zeros(K)];
    B2 = [-Dplus; zeros(K)];
    C =  [eye(N); zeros(K, N)];
    
    out = A1*x +  A2*xh + B1*ub(x) + B2*ub(xh) + C*w;
end 

function out = decomp2(x, xh, wh)
    % Backard-time decompisition function for backup dynamics
    global beta N K D
    
    Dplus = (D > 0).*D;  % Get positive part of D
    Dminus = (D < 0).*D; % Get negitive part of D

    A1 = [beta.*eye(N), zeros(N, K); Dminus', zeros(K)];
    A2 = [zeros(N), zeros(N, K); Dplus', zeros(K)];
    B1 = [Dplus; zeros(K)];
    B2 = [Dminus; zeros(K)];
    C =  [eye(N); zeros(K, N)];
    
    out = - A1*x - A2*xh + B1*ub(x) + B2*ub(xh) - C*wh;
end 

function out = Eb(x)
    % Embedding dynamics for CL system under backup controller
    % output = E(x, xhat)  
    global W N K

    out = [decomp(x(1:N+K), x(N+K+1:2*(N+K)), W(:, 1)); ...
           decomp(x(N+K+1:2*(N+K)), x(1:N+K), W(:, 2))];
end

function out = Eb2(x)
    % Embedding dynamics for backward-time CL system under backup controller
    % output = E(x, xhat)  
    global W N K

    out = [decomp2(x(1:N+K), x(N+K+1:2*(N+K)), W(:, 1)); ...
           decomp2(x(N+K+1:2*(N+K)), x(1:N+K), W(:, 2))];
end

function [Tout, Phi] = PhiE2(T, x0, xh0)
    dt = .001;
    % Simulate embedding trajectory
    Phi = [x0; xh0];
    Tout = 0:dt:T;
    sz = size(Tout, 2);
    for t = 1:sz - 1
        x_next = Phi(:, t) + dt.*Eb2(Phi(:, t));
        Phi(:, t+1) = x_next;
    end
end

function [Tout, Phi] = PhiE(T, x0, xh0)
    global dt
    % Simulate embedding trajectory
    Phi = [x0; xh0];
    Tout = 0:dt:T;
    sz = size(Tout, 2);
    for t = 1:sz - 1
        x_next = Phi(:, t) + dt.*Eb(Phi(:, t));
        Phi(:, t+1) = x_next;
    end
end

function out = DPhiE(T, x0)
    global N K
    
    dx = .0002;
    [~, Phi] = PhiE(T, x0, x0);
    Phi_T =  Phi(:, end);
   
    out = zeros(2*(N+K), N+K);
    for i = 1:N+K
        q = zeros(N+K, 1);  q(i) = dx;
        [~, Phi] = PhiE(T, x0 + q, x0 + q);
        out(:, i) =   (1/dx)*(Phi(:, end) - Phi_T);
    end
end

% -------------
%  safe set and barrier
% -------------

function [P, Delta] = get_Safe_set()
    % This code numerically checks whether the given quadratic Lyapunov 
    % function is valid for the 3-cart example subject to noise by randomly 
    % sampling on its boundary.
    fprintf('Getting Safe Set ... \n');
    global beta D Springs

    c = Springs(1)*Springs(2);
    k = -beta;
    springConst= Springs(1);
    stretch    = Springs(2);
    
    P = [c*eye(3) + D*D', k*D;...
         k*D',(k^2+c^2)*eye(2)+c*D'*D];

    FLAG = 0;
    Delta  = 1;
    while FLAG ~= 3
        Pinvsq = sqrtm(inv(P));
        n   = 2e5;
        r   = randn(5,n);
        r   = bsxfun(@rdivide,r,sqrt(sum(r.^2,1)));
        wthresh = .1*sqrt(3); %a 2-norm bound when the infinity-norm bound on the 3-dimensional disturbance is 0.075.
        Vdotcheck_extreme = -inf;

        for i=1:n
            x=r(:,i);
            y=sqrt(Delta)*Pinvsq*x;

            v=y(1:3); z=y(4:5);

            h=[springConst*tanh(stretch*z(1)); springConst*tanh(stretch*z(2))];
            ubackup = -D*h;

            vdot_nodist = -k*v+ubackup;
            zdot = D'*v;
            Vdot_temp = 2*y'*P;
            Vdotcheck = 2*y'*P*[vdot_nodist;zdot] + norm(Vdot_temp(1:3))*wthresh;

            Vdotcheck_extreme = max([Vdotcheck_extreme, Vdotcheck]);
        end
        Vdotcheck_extreme;
        if Vdotcheck_extreme > 0
            Delta = Delta + .25;
            FLAG = 0;
        else
            FLAG = FLAG + 1;
        end
    end
    fprintf('Safe Set verified for \n');
    Delta
end

function [S] = over_approx_Sb()
    fprintf('Approximating Safe Set with Rectangle... \n');
    global P L
    
    for i = 1:5
        cvx_begin quiet
            %Optimize over u
            variable x(5, 1)

            % QP Constraints
            L >= x' * P * x;

            minimize(x(i));
        cvx_end;
        S(i, 1) = x(i, 1);
    end
    for i = 1:5
        cvx_begin quiet
            %Optimize over u
            variable x(5, 1)

            % QP Constraints
            L >= x' * P * x;

            maximize(x(i));
        cvx_end;
        S(i, 2) = x(i, 1);
    end
end

function out = hsam(x)
    % This was numerically "verified" to be a robust Lyapunov function for
    % this particular L below.
    global  P L 
    
    % takes in full state - output is evaluation of barrier
    out = L - x.'*P*x;
end

function out = Dhsam(x)
global  P 
    % takes in full state - output is evaluation of barrier
    out=-2*x.'*P;
end

function out = hsamsoftmin(x)
    %returns the softmin of the quadratic barrier hsam evaluated on all
    %corners defined by the embedding state x
    m = 50; %also appears below
    [l, ~] = size(x);
    xlow   = x(1:l/2);
    xhigh  = x((1+l/2):end);
    corners= get_corners(xlow,xhigh);
    hs     = cellfun(@hsam, num2cell(corners,1));
    out    = -1/m*logsumexp(-m*hs);
end

function out = hsamsoftmingrad(x)
    m=50; %also appears above
    [l,~]=size(x);
    xlow=x(1:l/2);
    xhigh=x((1+l/2):end);
    corners=get_corners(xlow,xhigh);
    Dcorners=get_Dcorners(xlow,xhigh);
    hs=cellfun(@hsam,num2cell(corners,1));
    n2=length(hs);
    os=max(-m*hs);
    denominator=sum(exp(-m*hs-os));
    out=zeros(l,1)';
    for j=1:n2
        out=out+exp(-m*hs(j)-os)*Dhsam(corners(:,j))*Dcorners{j};
    end
    out=1/denominator*out;
end

% -------------
%  other
% -------------

function out = get_corners(xlow,xhigh)
    %given a lower and upper corner point, this returns all corners (as a
    %matrix)
    [n, ~]=size(xlow);
    xblock=[xlow, xhigh];
    i = de2bi((0: 2^n-1)')+1;
    [n2, m2] = size(i);
    out=zeros(n,n2);
    for j=1:n2
        for k=1:n
            out(k,j)=xblock(k,i(j,k));
        end
    end
end

function out = get_Dcorners(xlow, xhigh)
    % returns derivative of the corners with respect to xlow and xhigh as a
    % cell array. We actually only use xlow and xhigh to get the size since
    % the derivatives of the corners are highly structured binary matrices.
    
    [n, ~]  = size(xlow);
    i       = de2bi((0: 2^n-1)');
    [n2, ~] = size(i);
    out     = cell(n2,1);
    for j=1:n2
        out{j}=[diag(1-i(j,:)), diag(i(j,:))];
    end
end

function out = plot_rect(x_u, x_o)
    rect = [x_u(1), x_o(1), x_o(1), x_u(1); ...
            x_u(2), x_u(2), x_o(2), x_o(2)];
    out(1) = patch(rect(1, :), rect(2, :), [1, 1, 1], 'EdgeColor', [1, 1, 1]);
    out(2) = patch(rect(1, :), rect(2, :), 'r', 'FaceAlpha', .2);
end

% FOR VIDEO
% get to correct subplot
function video_maker(in)
    global create_video myVideo
    if create_video == 1
        if in == 1
            %Set up video recorder
            myVideo = VideoWriter('video.mp4', 'MPEG-4');
            myVideo.FrameRate = 5;   % Default 30
            myVideo.Quality = 50;    % Default 75
            open(myVideo);
        elseif in == 2
            Frame = getframe(gcf);
            writeVideo(myVideo, Frame);
        elseif in == 3
            close(myVideo);
        end
    end
end

function ploot(in)
    if in == 1
        subplot(52, 2, 1:20)
    elseif in == 2
        subplot(52, 2, 32:72)
    elseif in == 3
        subplot(52, 2, 85:2:103)
    elseif in == 4
        subplot(52, 2, 86:2:104)
    end
end

function CART = MakeCart(scale)
    CART = [];
    cartlegnth = scale;
    wheel_rad = scale/5;

    x_off = - cartlegnth/3;
    for i = 0:pi/4:2*pi
        % wheel 1
        point = [x_off + wheel_rad*sin(i); wheel_rad*cos(i)];
        CART = [CART, point];
    end
    x_off = cartlegnth/3;
    for i = 0:pi/4:2*pi
        % wheel 2
        point = [x_off + wheel_rad*sin(i); wheel_rad*cos(i)];
        CART = [CART, point];
    end
    point1 = CART(:, end) + [2*wheel_rad; 0];
    point2 = point1 + [0; 4*wheel_rad];
    point3 = point2 + [-(2/3)*cartlegnth - 4*wheel_rad; 0];
    point4 = point3 + [0; -4*wheel_rad];
    point5 = point4 + [2*wheel_rad; 0];
    CART = [CART, point1, point2, point3, point4, point5];
end




