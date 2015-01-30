function [v, xtraj,utraj,info ] = trajopt_hybrid_juggle( )
w = warning('off','Drake:RigidBodyManipulator:ReplacedCylinder');
load('bp-mode-2.mat'); %get u,x,and bendyplant from the seed traj opt
e = .6;%.857; %coefficient of restitution
bj = BendyJugglerPlant(BendyPlant(l,mass,extra_mass,k,c,N),e);

v = bj.p.constructVisualizer();
%v.xlim = [-.05,l+.05];
v.xlim = [-2*l,2*l];
%v.ylim = [-.1,.15];
v.ylim = [-2*l,2*l];

%get these automatically from xtraj, utraj....
x0 = Point(bj.getStateFrame());
x0(1) = 1; %mode
traj0 = xtraj.eval(0);
x0(2:bj.p.N+1) = traj0(1:bj.p.N);
x0(bj.p.N+2) = .3*bj.p.l; %.298*bj.p.l;
x0(bj.p.N+3) = 10.*bj.ball_radius; %3.*bj.ball_radius; %calculate zero crossing time height?

x0(bj.p.N+5:2*bj.p.N+4) = traj0(bj.p.N+1:end);
%x0(2*bj.p.N+6) = -1.; %q1dot

%u0 = utraj;
%u0 = setOutputFrame(u0,bj.p.getInputFrame);
%[xtraj,utraj] = simulate(cascade(u0,bj),[0;1.5],x0);
%v.draw(0,xtraj.eval(0));
%v.draw(0,xtraj.eval(.105));
%v.draw(0,xtraj.eval(.219));

%info=0;
%playback(v,xtraj,struct('slider',1)); 

%return;

N = [3;3]; %time steps for each mode
mode_sequence = [1;2]; %to enforce that a guard is tripped
durations = {[.02 .5],[.02 .5]}; %upper and lower bounds for duration of each mode in mode_sequence
options.u_const_across_transitions = true;
options.periodic = true;
traj_opt = HybridTrajectoryOptimization(@DircolTrajectoryOptimization,bj,mode_sequence,N,durations,options);

%extracted from seed trajectory
x0 = x0(1:end-1);
t1 = .256;
x1 = [-1.20766497257925;0.0318361731825689;0.144529970484959;0.0956688506813374;-0.0219308384812992;0.0541011215893433;0.0991117898432538;-0.0210902882975062;-0.021792154320133;0.0756391267046635;0.250003387345963;0.109165444990245;-0.0210043670538273;-0.0317376826557072;0.458885699429571;-0.598335958701967;0.137662463804818;0.280635399528564;-0.312683113949808;-0.1135260074932;0.273320399590499;0.0240788401921236;-0.22143644185827;-5.70781506476142e-05;0.54822703149028;4.15767959978957e-07];
tf = .52;
xf = x0;
t_init{1} = linspace(0,t1,N(1));
t_init{2} = linspace(0,tf-t1,N(2));

traj_init{1}.u = PPTrajectory(foh(t_init{1},utraj.eval(t_init{1})));
traj_init{1}.u = traj_init{1}.u.setOutputFrame(bj.p.getInputFrame);
traj_init{2}.u = PPTrajectory(foh(t_init{2},utraj.eval(t_init{2})));
traj_init{2}.u = traj_init{2}.u.setOutputFrame(bj.p.getInputFrame);
traj_init{1}.x0 = x0;
traj_init{2}.x0 = x1;

%n_var = size(x0,1);
%traj_opt = traj_opt.addModeStateConstraint(1,BoundingBoxConstraint(-inf(n_var,1),inf(n_var,1)),1);
traj_opt = traj_opt.addModeStateConstraint(1,BoundingBoxConstraint(.5*x0(bj.p.N+3),2*x0(bj.p.N+3)),1,bj.p.N+2);
%traj_opt = traj_opt.addModeStateConstraint(1,ConstantConstraint(x0(bj.p.N+2)),1,bj.p.N+1);

%traj_opt = traj_opt.addModeStateConstraint(2,BoundingBoxConstraint(x0(bj.p.N+3),inf),N(1),bj.p.N+3);

traj_opt = traj_opt.addModeRunningCost(1,@cost);
traj_opt = traj_opt.addModeRunningCost(2,@cost); %HACK!

traj_opt = traj_opt.compile();
traj_opt = traj_opt.setCheckGrad(true);
%snprint('snopt.out'); 
%disp(traj_opt);
tic
[xtraj,utraj,z,F,info] = solveTraj(traj_opt,t_init,traj_init);
toc
playback(v,xtraj,struct('slider',1)); 
end

function [g,dg] = cost(~,x,u)
    n = size(x,1)/2; %num positions
    Q = diag([0.,0.*ones(1,n-1), 0., 0.*ones(1,n-1)]);     %state cost
    R = 0.0001;                                             %effort cost
    g = x'*Q*x + u'*R*u;
    dg = [0, 2*x'*Q, 2*u'*R];
    %disp(size(g));
    %disp(size(dg));
end

% function [g,dg] = final_cost(T,x)
%     n = size(x,1)/2;
%     Q = diag([0.,zeros(1,n-1), -2., zeros(1,n-1)]);     %state cost
%     g = x'*Q*x;
%     dg = [0, 2*x'*Q];
% end