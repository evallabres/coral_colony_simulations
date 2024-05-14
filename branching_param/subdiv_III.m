
   %      o           o
   %     / \         / \
   %    x   x  ---> o---o
   %   /     \     / \ / \
   %  o---x---o   o---o---o
   
 function [U14,F14,E14,B14] = subdiv_III(offset,V,F,B,N)
     %fprintf('subdivIII')
     E14 = [F(:,2) F(:,3);F(:,3) F(:,1);F(:,1) F(:,2)];
     U14 = smooth_midpoint(V(E14(:,1),:), V(E14(:,2),:), N(E14(:,1),:), N(E14(:,2),:));
     % indices of midpoints
     nu = size(U14,1);
     i1 = offset+(1:(nu/3))';
     i2 = offset+((nu/3)) + (1:(nu/3))';
     i3 = offset+((nu/3)+(nu/3)) + (1:(nu/3))';
     % new face indices, 4 new faces for each original face.
     F14 = [F(:,1) i3 i2 ; F(:,2) i1 i3 ; F(:,3) i2 i1 ; i1 i2 i3];
     B14 = [B(:);B(:);B(:);B(:)];
   end