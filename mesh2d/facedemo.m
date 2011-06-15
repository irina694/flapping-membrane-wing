function facedemo(n)

% FACETEST: Example polygonal geometries for MESHFACES.
%
%  facedemo(n)
%
% Darren Engwirda - 2007

switch n
   case 1

      close all
      %node = [0.0, 0.0; 1.0,0.0; 1.0,1.0; 0.0,1.0; 1.01,0.0; 1.01,1.0; 3.0,0.0; 3.0,1.0];
      node = [0.0, 0.0; 1.0,0.0; 1.0,1.0; 0.0,1.0; 2,0.0; 3.0,0.0; 3.0,1.0];
      %edge = [1,2; 2,3; 3,4; 4,1; 2,5; 5,6; 6,3; 5,7; 7,8; 8,6];
      edge = [1,2; 2,3; 3,4; 4,1; 2,5; 5,3; 5,6; 6,7; 7,3];
      face{1} = [1,2,3,4];
      face{2} = [5,2,6];
      face{3} = [8,7,6,9];

      hdata.hmax  = .2;
      options.maxit = 20;
      options.dhmax = 0.2;
      
      [p,t,fnum] = meshfaces(node,edge,face,hdata,options);

      hold on
      i = 3;
      foo = cell2mat(face(i));
      for j = 1:length(foo)
          plot([node(edge(foo(j),1),1),node(edge(foo(j),2),1)],[node(edge(foo(j),1),2),node(edge(foo(j),2),2)],'k.')
      end
      
   case 2

      % Geometry
      dtheta = pi/36;
      theta = (-pi:dtheta:(pi-dtheta))';
      node1 = [cos(theta), sin(theta)];
      node2 = [-2.0,-2.0; 2.0,-2.0; 2.0,2.0; -2.0, 2.0];
      edge1 = [(1:size(node1,1))',[(2:size(node1,1))'; 1]];
      edge2 = [1,2; 2,3; 3,4; 4,1];

      edge = [edge1; edge2+size(node1,1)];
      node = [node1; node2];

      face{1} = 1:size(edge1,1);
      face{2} = 1:size(edge,1);

      meshfaces(node,edge,face);

   otherwise
      error('Invalid demo. N must be between 1-2');
end

end      % facetest()
