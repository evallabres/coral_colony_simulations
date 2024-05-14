  
   %      o           o
   %     / \         /|\
   %    /   \  ---> / | \
   %   /     \     /  |  \
   %  o---x---o   o---o---o
   
   function [U12,F12,E12,B12] = subdiv_I(offset,V,M12,B,N)
   % M12 is a nx3 matrix similar to faces, but it should have only the faces
   % that i subdivide in two, and the two first indices of each row
   % correspond to the split edge. 
        %fprintf('subdivI')
        E12 = M12(:,1:2);
        U12 = smooth_midpoint(V(M12(:,1),:), V(M12(:,2),:), N(M12(:,1),:), N(M12(:,2),:));
        nu = size(U12,1);
        i1 = offset + (1:nu)';
        O1 = M12(:,1);
        O2 = M12(:,2);
        O3 = M12(:,3);
        F12 = [O3 O1 i1; i1 O2 O3];
        B12 = [B(:);B(:)];
   end

   
  