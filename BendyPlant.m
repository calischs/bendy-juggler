classdef BendyPlant < PlanarRigidBodyManipulator

    properties
        l;              %total length
        mass;           %total mass
        extra_mass;     %extra mass at end
        k;              %joint stiffness
        c;              %joint damping
        N;              %number of links
    end
    
  methods 
    function obj = BendyPlant(l,mass,extra_mass,k,c,N)
      commandStr = sprintf('python write_bendy_urdf.py %f %f %f %f %f %d',l,mass,extra_mass,k,c,N);
      [status, commandOut] = system(commandStr); 
      assert(status==0); %make sure system call to python completes
      obj = obj@PlanarRigidBodyManipulator('bendy.xml');
      obj.l=l; obj.mass=mass; obj.extra_mass=extra_mass; obj.k=k; obj.c=c; obj.N=N;
      %obj = setOutputFrame(obj,getStateFrame(obj)); %do we need this?
    end
    
    function [pts,ds0,ds1] = segment_endpoints(obj,x,i)
        %return array with columns as end points of segment i
        %ds0, ds1 have size (2,nq)
        assert(size(x,1)==2*obj.num_positions); %make sure we're not talking about the hybrid system.
        seg_num = i+1; %discount base frame!
        kinsol = obj.doKinematics(x(1:obj.num_positions));
        [pt1,ds0] = obj.forwardKin(kinsol,seg_num,[0.;0.]);
        [pt2,ds1] = obj.forwardKin(kinsol,seg_num,[0.;-obj.segment_length()]); %end point in link coords
        pts = [pt1 pt2];
    end
    
    function dl = segment_length(obj)
        %return length of segments
        dl = obj.l/obj.N;
    end
    function S = saveobj(obj,xtraj,utraj)
      % Save property values in struct
      % Return struct for save function to write to MAT-file
         S.l = obj.l;
         S.mass = obj.mass;
         S.extra_mass = obj.extra_mass;
         S.k = obj.k;
         S.c = obj.c;
         S.N = obj.N;
         S.xtraj = xtraj;
         S.utraj = utraj;
    end
    
%     function x0 = getInitialState(obj)
%         x0 = Point(obj.getStateFrame());
%         x0.joint1 = -pi/2;
%         x0.joint2 = -.01;
%         x0.base_x = .5*bj.p.segment_length();
%         x0.base_z = bj.ball_radius;
%         x0.base_z = 1.5*bj.ball_radius;
%     end
  end
end