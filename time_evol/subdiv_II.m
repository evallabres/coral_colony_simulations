  
   
   %      o           o
   %     / \         /|\
   %    x   \  ---> o | \
   %   /     \     / \|  \
   %  o---x---o   o---o---o
   
   function [U13,F13,E13,B13] = subdiv_II(offset,V,M13,B,N)
   % M12 is a nx3 matrix similar to F, but the two pairs of vertices (1,2)
   % and (2,3), are the ones that should be subdivided. 
        %fprintf('subdivII')
        E13 = [M13(:,1:2); M13(:,2:3)];
        U13a = smooth_midpoint(V(M13(:,1),:), V(M13(:,2),:), N(M13(:,1),:), N(M13(:,2),:));
        U13b =smooth_midpoint(V(M13(:,2),:), V(M13(:,3),:), N(M13(:,2),:), N(M13(:,3),:));
        U13 = [U13a; U13b];
        nu = size(U13,1);
        i1 = offset+(1:(nu/2))';
        i2 = offset+((nu/2)) + (1:(nu/2))';
        O1 = M13(:,1);
        O2 = M13(:,2);
        O3 = M13(:,3);
        %F13 = [O1 i1 i2; i1 O2 i2; O1 i2 O3];
        F13 = [O1 i1 O3; i1 i2 O3; i1 O2 i2];
        B13 = [B(:);B(:);B(:)];
   end
   
     