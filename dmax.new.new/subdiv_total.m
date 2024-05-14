
 function [VV,FF,BB] = subdiv_total(V,F,B,N,dmax)
 
     [E12,F12,V12,M12,B12,E13,F13,V13,M13,B13,E14,F14,V14,M14,B14,F11,B11,iaux,ht] =deal([]);
     nf = size(F,1);
     n = size(V,1);
     %%useful definitions
     iedge = [1 2; 3 1; 2 3;];
     icyc = [1 2 3];
     [iaux(3,:),iaux(4,:),iaux(5,:)] = deal([3 1 2],[1 2 3],[2 3 1]);

     for f = 1:nf
        [g,h] = deal([]);
        E = [F(f,iedge(1,:));F(f,iedge(2,:));F(f,iedge(3,:))];
        Vdif = V(E(:,1),:)-V(E(:,2),:);
        d = [norm(Vdif(1,:)) norm(Vdif(2,:)) norm(Vdif(3,:))];
        [isub, nsub] = deal(find(d > dmax), sum(d > dmax));
        if nsub == 3
            M14 = [M14; F(f,:)];
            B14 = [B14; B(f)];
        elseif nsub == 2
            naux = sum(isub);
            h = iaux(naux,:);  % this could be without the auxiliary variable.        
            M13 = [M13; F(f,h)];
            B13 = [B13; B(f)];
        elseif nsub == 1
            h = circshift(icyc, isub-1);
            M12 = [M12; F(f,h)];
            B12 = [B12; B(f)];             
        else 
            F11 = [F11; F(f,:)];
            B11 = [B11; B(f)];
        end        
     end
     %% FUSIONS
     %%% fusionI
     if isempty(ht) == 0 
        Vf1 = (V(ht(:,1),:)+V(ht(:,2),:))/2
         [V(ht(:,1),:), V(ht(:,2),:)]= deal([Vf1]);
     end
     %% SUBDIVISIONS
     if isempty(M14) == 0 
        [V14,F14,E14,B14] = subdiv_III(n,V,M14,B14,N);
     end
     if isempty(M13) == 0 
        [V13,F13,E13,B13] = subdiv_II(n+size(V14,1),V,M13,B13,N);
     end
     if isempty(M12) == 0 
        [V12,F12,E12,B12] = subdiv_I(n+size(V14,1)+size(V13,1),V,M12,B12,N);
     end
     
     %% FINAL GEOMETRY
     %% final faces and branches
     FF = [F14;F13;F12;F11];
     BB = [B14;B13;B12;B11];
     %% new vertices
     U = [V14;V13;V12];
     %% eliminating double vertices due to subdivision
     %%% find the edges you have split twice
     EU = [E14;E13;E12];
     [~,I,J] = unique(sort(EU,2),'rows');
     %%% relabel vertices and append old to new
     U = U(I,:);
     VV = [V ; U]; 
     %%% relabel faces
     J = [(1:n)';n+J];
     FF = J(FF);
     BB = BB;

 end
