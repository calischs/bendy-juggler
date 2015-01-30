classdef BendyJugglerPlant < HybridDrakeSystem
    properties
       p; %PlanarRigidBodyManipulator 
       ball_radius; %radius of the ball (must match urdf)
       e; %coefficient of restitution (in [0,1])
    end
    methods 
        function obj = BendyJugglerPlant(bp,e)
          obj = obj@HybridDrakeSystem(1,2*bp.N+6); %num_inputs, num_outputs
          %add ball
          bp = bp.addRobotFromURDF('ball.urdf',zeros(3,1),zeros(3,1),struct('floating',true));
          obj.p = bp; %keep PRBM around
          obj.ball_radius = .03; %radius of ball -- must match urdf.  TODO: better ball
          obj.e = e; %coefficient of restitution
          
          obj = setInputFrame(obj,obj.p.getInputFrame);
          obj = setOutputFrame(obj,obj.p.getOutputFrame);

          [obj,flight_mode_ind] = obj.addMode(obj.p,'flight'); %flight
          [obj,flight_mode2_ind] = obj.addMode(obj.p,'flight2'); %HACK, shouldn't really need another mode
          
          % logical guard minimum
          guard =  notGuard(obj, @(obj,t,x,u) obj.segment_distance_func(t,x,u,1,obj.ball_radius) );
          for i=2:bp.N
           guard = andGuards(obj,guard,notGuard(obj, @(obj,t,x,u) obj.segment_distance_func(t,x,u,i,obj.ball_radius) ));
          end
          obj = obj.addTransition(flight_mode_ind,notGuard(obj,guard),@collisionDynamics,false,true); 
          obj = obj.addTransition(flight_mode2_ind,notGuard(obj,guard),@collisionDynamics,false,true); %HACK, shouldn't really need another mode
          
          % smooth exp minumum approx
          %for i=1:bp.N
          %   guards{i} = @(obj,t,x,u) obj.segment_distance_func(t,x,u,i,obj.ball_radius);
          %end
          %guard = obj.smooth_min_exp(guards,1e2);
          %obj = obj.addTransition(flight_mode_ind,guard,@collisionDynamics,false,true); 
          %obj = obj.addTransition(flight_mode2_ind,guard,@collisionDynamics,false,true); %HACK, shouldn't really need another mode
                      
          obj = setSimulinkParam(obj,'InitialStep','1e-4','MaxStep','1e-1');
        end
        
        function [pt,dptdq] = ball_pos(obj,xx)
            np= obj.p.N;
            kinsol = obj.p.doKinematics(xx(1:obj.p.num_positions));
            [pt,dptdq] = obj.p.forwardKin(kinsol,np+4,[0.;0.]);
            %pt = [xx(np+1); xx(np+2)];
            %if nargout>1
            %   I = eye(np+3);
            %   dptdq = I([np+1,np+2],:); 
            %end
        end
    
        function [xp,mode,status,dxp] = collisionDynamics(obj,mode,t,xm,u)
          %disp('collision dynamics called');
          %disp(size(u));
          if (mode==1) mode=2;  % switch modes, HACK
          else mode=1; end   
            
          np = obj.p.num_positions;  %in whole system
          qm = xm(1:np); vm = xm(np+1:2*np);
          %get inertia matrix H and its derivative wrt x
          
          [H,~,~,dHdx,~,~] = obj.p.manipulatorDynamics(qm,vm);
          dHdqm = dHdx(:,1:np); %discard dHdvm=0;
          
          %get kinsol for specific configuration
          %kinsol = obj.p.doKinematics(qm);   
          kinsol_dumb = obj.p.doKinematics(qm,true,false); %non-mex version  

          %get contact jacobian on manipulator
          [i,pt_on_manip,direction,dpt_on_manipdq,ddirectiondq] = obj.closest_point_on_manip_to_ball(t,xm,u);
          [~,am,damdqm] = obj.p.forwardKin(kinsol_dumb,i,pt_on_manip,0);
          %adjust damdqm for dpt_on_manipdq
          T23 = obj.p.T_2D_to_3D; T32 = T23';
          dTdq = kinsol_dumb.dTdq{i}; %adapted from forwardKin
          cross_term = zeros(2,np^2);
          for j=1:np
              cross_term(:,np*(j-1)+1:np*j) = T32*dTdq(3*(j-1)+1:3*j,1:3)*T23*dpt_on_manipdq;
          end
          damdqm = damdqm + cross_term;
          damdqm = reshape(damdqm'*direction,np,np) + ddirectiondq'*am;
          am = am'*direction; %pick out the jacobian in the normal direction.
          
          %get contact jacobian on ball
          ball_body_num = obj.p.N+3; %13; %BAD
          pt_on_ball = [0;0]; %this should be at the contact point, but w/o friction, I don't think it matters...
          [~,ab,dabdqm] = obj.p.forwardKin(kinsol_dumb,ball_body_num,pt_on_ball,0);
          dabdqm = reshape(dabdqm'*direction,np,np) + ab'*ddirectiondq;
          ab = ab'*direction;

          a = ab-am; dadqm = dabdqm-damdqm;
          
          Hxa = inv(H)*a;%H\a; %save inversion of matrix.
          Hnorm = a'*Hxa;
          Hi = inv(H);
          L = -(1+obj.e)*(a'*vm)/Hnorm; %collision impulse
          vp = vm + Hxa*L; %velocity transition
          %new state has same position, new velocity
          
          dt = 0;%1e-5; %add a small dt*vp to the position.
          xp = [xm(1:np)+dt*vp; vp]; 

          if (nargout>3)
              %do dxpdq
              %maybe should include the small dt in here too?
              dxp = zeros(2*np,2*np);
              dxp(1:np,1:np) = eye(np,np); %dqpdqm
              M = eye(np,np)-(1+obj.e)*(Hxa*(a'))/Hnorm;
              
              dMdqm_vm = zeros(np,np);
              dHdqm_chunksize = size(dHdqm,2); %kind of a hack
              for j=1:np
                d1 = (-1/Hnorm^2)*(2*Hxa'*dadqm(:,j) - Hxa'*dHdqm(dHdqm_chunksize*(j-1)+1:dHdqm_chunksize*j,:)*Hxa)*Hxa*a'*vm;
                d2 = -(1/Hnorm)*Hi*dHdqm(dHdqm_chunksize*(j-1)+1:dHdqm_chunksize*j,:)*Hxa*a'*vm;
                d3 = (1/Hnorm)*Hi*(dadqm(:,j)*a'*vm + a*dadqm(:,j)'*vm);
                dMdqm_vm(:,j) = -(1+obj.e)*(d1 + d2 + d3);
              end
              dxp(np+1:2*np, np+1:2*np) = M; %dvpdvm
              dxp(np+1:2*np, 1:np) = dMdqm_vm; %dvpdqm
              
              dxp = [zeros(2*np,2), dxp, zeros(2*np,1)]; %don't forget dxpdt and dxpdu, copy CompassGaitPlant, not sure on 2
          end
          status = 0;
        end
        
        function [dist,d_dist] = segment_distance_func(obj,~,x,~,i,radius)
        %return distance from ball center to segment i
            if nargout>1
                [s,ds0dq,ds1dq] = obj.p.segment_endpoints(x,i);
                [p,dpdq] = obj.ball_pos(x);
                [dist,ddds0,ddds1,dddp] = obj.seg_dist( s, p );
                ddistdq = dddp*dpdq + ddds0*ds0dq + ddds1*ds1dq; %chain rule
                d_dist = [0, ddistdq, zeros(size(ddistdq)), 0]; %add placeholders for t,v, and u
            else
                dist = obj.seg_dist( obj.p.segment_endpoints(x,i), obj.ball_pos(x) );
            end
            dist = dist - radius^2;
            %dist = sqrt(dist) - radius; %smooth exp min likes actual
            %distance, instead of squared...need to take this into account
            %in derivatives...
        end
        
        function [i,proj,dir,dprojdq,ddirdq] = closest_point_on_manip_to_ball(obj,t,x,u)
            %find segment i with closest point to ball
            %find point proj in segment i coordinates of this point
            %find direction dir from ball center to closest point on
            %manipulator
            [pt,dpdq] = obj.ball_pos(x);
            dists = zeros(obj.p.N,1);
            for j=1:obj.p.N
                dists(j,1) = obj.segment_distance_func(t,x,u,j,0);
            end
            [~,i] = min(dists,[],1); %get argmin
            
            [s,ds0dq,ds1dq] = obj.p.segment_endpoints(x,i);
            [alpha,dir,dalphads0,dalphads1,dalphadp,ddirds0,ddirds1,ddirdp] = obj.seg_closest_point(s,pt);
            ddirdq = ddirdp*dpdq + ddirds0*ds0dq + ddirds1*ds1dq;
            dalphadq = dalphadp*dpdq + dalphads0*ds0dq + dalphads1*ds1dq;
            l = obj.p.segment_length();
            proj = [0; -alpha*l]; %in original segment coords
            dprojdq = [zeros(size(dalphadq)); -l*dalphadq];
        end
        
        
        function [alpha,dir,dalphads0,dalphads1,dalphadp,ddirds0,ddirds1,ddirdp] = seg_closest_point(obj,s,p)
                %compute parameter alpha along the segment that is closest to pt
                %0 -> first end point, 1-> second end point
                %also compute direction from pt to closest point
                %s has columns with segment end points
                U = p-s(:,1); V = s(:,2)-s(:,1);
                uv = U'*V; vv = V'*V;
                if uv<0
                    alpha = 0.; %s(:,1) is closest
                    dalphadp = zeros(1,2);
                    dalphads0 = zeros(1,2);
                    dalphads1 = zeros(1,2);
                    dir = U;
                    ddirdp = eye(2);
                    ddirds0 = -eye(2);
                    ddirds1 = zeros(2);
                elseif uv>=vv
                    alpha = 1.; %s(:,2) is closest
                    dalphadp = zeros(1,2);
                    dalphads0 = zeros(1,2);
                    dalphads1 = zeros(1,2);
                    dir = p-s(:,2);
                    ddirdp = eye(2);
                    ddirds0 = zeros(2);
                    ddirds1 = -eye(2);
                else
                    alpha = uv/vv;
                    dalphadp = (1/vv)*V';
                    dalphads0 = 2*(uv/vv^2)*V' - (1/vv)*(U'+V');
                    dalphads1 = (1/vv)*U' - 2*(uv/vv^2)*V';
                    dir = U-alpha*V;
                    ddirdp = eye(2) - V*dalphadp;
                    ddirds0 = -eye(2) + alpha*eye(2) - V*dalphads0;
                    ddirds1 = -alpha*eye(2) - V*dalphads1;
                end
                %normalize direction and adjust derivatives
                %[dist, ddds0,ddds1,dddp] = obj.seg_dist(s,p);
                %ddirdp = (dist*ddirdp - dir*dddp)/dist^2;
                %ddirds0 = (dist*ddirds0 - dir*ddds0)/dist^2;
                %ddirds1 = (dist*ddirds1 - dir*ddds1)/dist^2;
                %dir = dir/dist;
        end
    end
             
    methods (Static=true)
        
        function [d,ddds0,ddds1,dddp] = seg_dist(s,p)
            %s has columns with segment end points
            U = p-s(:,1); V = s(:,2)-s(:,1);
            uv = U'*V; vv = V'*V; uvovv = uv/vv;
            if uv<0
                d = U'*U;
            elseif uv>=vv
                d = (p-s(:,2))'*(p-s(:,2));
            else
                d = (U-uvovv*V)'*(U-uvovv*V);
            end
            if nargout>1
                if uv<0
                    dddp = 2*U';
                    ddds0= -2*U';
                    ddds1= zeros(1,2);
                elseif uv>=vv
                    dddp = 2*(p-s(:,2))';
                    ddds0= zeros(1,2);
                    ddds1= -2*(p-s(:,2))';
                else
                    dddp = 2*(U-uvovv*V)'*(eye(2)-1/vv*(V*V'));
                    ddds0= 2*(U-uvovv*V)'*((uvovv-1)*eye(2) + 1/vv*(V*U'+V*V') - 2*uvovv/vv*(V*V'));
                    ddds1= 2*(U-uvovv*V)'*(-uvovv*eye(2) - 1/vv*(V*U') + 2*uvovv/vv*(V*V'));
                end    
            end
        end
        
        
        
        function g = smooth_min_exp(guards,k)
            %smooth minimum combination of guards: http://www.iquilezles.org/www/articles/smin/smin.htm
            %guards is a 1xn cell array of functions for each guard
            %k controls smoothness vs. accuracy.  Large k is smoother
            n = size(guards,2);
            function [z,dz] = new_guard(obj,t,x,u)
                gi = zeros(n,1);
                if nargout>1
                    dgi = zeros(n,1+size(x,1)+1); %remember ddt, ddu
                end
                for i=1:n
                    guard = guards{i};
                    if nargout==1
                        [gi(i,1),~] = guard(obj,t,x,u);
                    else
                        [gi(i,1),dgi(i,:)] = guard(obj,t,x,u);
                    end
                end
                
                z = -log(sum(exp(-k*gi),1))/k; 
                if nargout>1
                    dz = sum(exp(-k*gi)'*dgi,1) ./ sum(exp(-k*gi),1);
                end
            end
            g = @new_guard; 
        end
    end
end