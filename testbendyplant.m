function [v, xtraj,utraj,info ] = testbendyplant( )

w = warning('off','Drake:RigidBodyManipulator:ReplacedCylinder');
bp = BendyPlant(.5,.2,1.,.5e2,0.1,10); %(l,mass,extra_mass,k,c,N)
v = bp.constructVisualizer();
%v.xlim = [-.05,.55];
v.xlim = [-1,1.];
%v.ylim = [-.15,.15];
v.ylim = [-.5,.5];
x0 = zeros(2*bp.N,1);
x0(1) = -pi/2;

%find holding torque to make a straight-out fixed point
%figure(1);
%v.draw(0,x0);
holding_torque = -(bp.l*9.81*bp.extra_mass + .5*bp.l*9.81*bp.mass);
%disp(holding_torque);
fpp = FixedPointProgram(bp);
fpp = fpp.addConstraint( ConstantConstraint(-pi/2),1);
[x_star,u_star] = fpp.findFixedPoint(x0,holding_torque);
holding_torque= u_star.double;

%v.draw(0,x_star)
%ylim([-.15,.05])
%disp(holding_torque);
%return;
%construct control trajectory manually
%u0 = FunctionHandleTrajectory(@(t) holding_torque,1,[0 5]);
%u0 = setOutputFrame(u0,bp.getInputFrame);
%[xtraj,utraj] = simulate(cascade(u0,bp),[0;5.],x0);
%v.playback(xtraj,struct('slider',true));

%construct the same control trajectory with a Constant frequency trajectory
%amp  = 10.;
%controller = SinusoidalController(bp,holding_torque,amp);
%freq_sys = cascade(controller,bp);
%w =  FunctionHandleTrajectory(@(t) .1, 1,[0 5]);
%w = setOutputFrame(w,freq_sys.getInputFrame);
%xtraj = simulate(cascade(w,freq_sys),[0;5.],x0);

%construct the same control trajectory with a Constant frequency trajectory
%amp  = 10.;
%fdbp = FrequencyDrivenBendyPlant(bp,holding_torque,amp);
%fdbp = FrequencyDrivenBendyPlant(bp,0.,amp);
%w =  FunctionHandleTrajectory(@(t) .1, 1,[0 5]);
%w = setOutputFrame(w,fdbp.getInputFrame);
%x0(1) = 0;
%xtraj = simulate(cascade(w,fdbp),[0;5.],x0);


traj_opt = true;
if traj_opt

    %freq_sys = FrequencyDrivenBendyPlant(bp,holding_torque,10.);
    %freq_sys = FrequencyDrivenBendyPlant(bp,0.,10.);

    N = 20; %n time points in traj opt
    prog = DircolTrajectoryOptimization(bp,N,[.1,2]);
    % all q at t = 0 are 0
    prog = prog.addStateConstraint(ConstantConstraint(-pi/2),1,1); %constraint, time index, x_indices
    prog = prog.addStateConstraint(ConstantConstraint(zeros(1,bp.N-1)),1,2:bp.N);

    %q1 at t=N/2 is a specific value
    %prog = prog.addStateConstraint(ConstantConstraint(-pi/2+0.05),N/2,1);
    
    % all q and qdot at t = 0 equal all qdot at t = N
    ns = bp.getNumStates;
    A = [diag([ones(ns/2,1); ones(ns/2,1)]), diag([-ones(ns/2,1); -ones(ns/2,1)])];
    periodic_constraint = LinearConstraint(zeros(ns,1),zeros(ns,1),A);
    prog = prog.addStateConstraint(periodic_constraint,{[1,N]});
    
    % u at t = 0 equals u at t = N
    prog = prog.addInputConstraint(LinearConstraint(zeros(2,1),zeros(2,1),eye(2)),{[1,N]});
  
    
    prog = prog.addRunningCost(@cost);
    
    %q1dot at t=0 is greater than a specific value
    %A = zeros(1,ns); A(ns/2+1)=1.;
    %prog = prog.addStateConstraint(LinearConstraint(.5,inf,A),{1});
    %prog = prog.addFinalCost(@final_cost);
    
    %prog = prog.addTrajectoryDisplayFunction(@(t,x,u)plotDircolTraj(t,x,u));
    
    
    ht = ConstantTrajectory(holding_torque);
    ht = setOutputFrame(ht,bp.getInputFrame);
    initial_guess.u = ht;
    initial_guess.x = ConstantTrajectory([-pi/2;zeros(2*bp.N-1,1)]);
    ts = linspace(0,2.,20);
    omega_0 = 1.*2*pi;
    %initial_guess.x = PPTrajectory(foh(ts,...
    %[-pi/2+.1*sin(omega_0*ts); -.01*ones(bp.N-1,1)*sin(omega_0*ts);...
    %.1*omega_0*cos(omega_0*ts); -.01*omega_0*ones(bp.N-1,1)*cos(omega_0*ts)]));
    prog = prog.setCheckGrad(true);

    [xtraj,utraj] = prog.solveTraj(ts(end),initial_guess);
    
    %plot_traj(utraj);
    %plot_states(xtraj);
    
    bpsave = bp.saveobj(xtraj,utraj);
    disp(bpsave);
    save('bp-mode-2.mat','-struct','bpsave');
    
    t = linspace(xtraj.tspan(1),xtraj.tspan(2),500);
    figure(24);
    x = xtraj(1).eval(t);
    plot(x(1,:),x(bp.N+1,:)); 
    %ylim([-1,1]);
    title('q_0 phase portrait');
    
    %play a couple cycles
    T = xtraj(1).tspan(2);
    looptraj = xtraj;
    for i=1:4
        looptraj = looptraj.append(xtraj.shiftTime(i*T));
    end
    %v.draw(0,xtraj.eval(0));
    %v.draw(0,xtraj.eval(.861));
    v.playback(looptraj,struct('slider',true));
    beep;
end
end

function [g,dg] = cost(dt,x,u)
    n = size(x,1)/2;
    Q = diag([1.,0.*ones(1,n-1), 0., 1.0*ones(1,n-1)]);     %state cost
    R = .0;                                             %effort cost
    g = x'*Q*x + u'*R*u;
    dg = [0, 2*x'*Q, 2*R*u];
    
    %add term minimizing the end effector's z displacement from zero
    %[G,dG] = z_disp(x);
    %g = g + G;
    %dg = dg + dG;
end

function [g,dg] = final_cost(T,x)
    n = size(x,1)/2;
    Q = zeros(2*n,2*n);%diag([0.,zeros(1,n-1), -2., zeros(1,n-1)]);     %state cost
    g = x'*Q*x + T;
    dg = [1, 2*x'*Q];
end

function [g,dg] = z_disp(x)
    %z displacement squared of end effector, and derivative wrt (t,x,u)
    n = size(x,1)/2;
    l = 10*.5/n;                   %link length, stupid
    num_back = 0;               %number of joints to count back from end
    ths = cumsum(x(1:n));
    displace = sum(l*cos(ths(1:(end-num_back))),1); %experiment to get second to last joint
    g = displace.*displace;
    if nargout>1
        sins = cumsum(-l*sin(ths(end:-1:1)));
        d_displace = [sins(end:-1:(1+num_back))',zeros(1,num_back), zeros(1,n)]; %remember qdots
        dg = [0, 2*displace*d_displace, 0]; %remember t and u
    end
end

function plotDircolTraj(t,x,u)
  figure(25);
  n = size(x,1);
  h=plot(x(1,:),x(n/2+1,:),'r.-');
  drawnow;
  delete(h);
end

function plot_traj(traj)
    t = linspace(traj.tspan(1),traj.tspan(2),50);
    figure(24);
    plot(t,traj.eval(t));  
    title('Control trajectory');
end
function plot_states(traj)
    t = linspace(traj(1).tspan(1),traj(1).tspan(2),50);
    figure(23);
    x = traj.eval(t);
    hold on;
    plot(t,x(1,:)',t,x(2,:)',t,x(3,:)',t,x(4,:)',t,x(5,:)',x(6,:)',x(7,:)',x(8,:)',x(9,:)','Color','b');
    plot(t,x(9+1,:)',t,x(9+2,:)',t,x(9+3,:)',t,x(9+4,:)',t,x(9+5,:)',x(9+6,:)',x(9+7,:)',x(9+8,:)',x(9+9,:)','Color','g');
    hold off;
    xlim([t(1),t(end)]);
    title('State trajectories');
    figure(22);
    zs = zeros(size(t));
    for i=1:size(t,2)
        zs(i) = z_disp(traj.eval(t(i)));
    end
    plot(t,zs);  
    title('z displacements');
end