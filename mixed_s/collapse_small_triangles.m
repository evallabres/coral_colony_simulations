
%to make the function faster, we should out only the label of the faces
%that are not eliminated, and then we can translate be in the main code.

function [FF,indexv] = collapse_small_triangles(V, F, dblA, min_dblarea)
    iedge = [1 2; 2 3; 3 1;];
    
    %bbd = bounding_box_diagonal(V);
    % Minimum area tolerance
    %min_dblarea = 2.0 * eps * bbd^2;
    
    FIM = (1:size(V, 1))';
    num_edge_collapses = 0;
    
    % Loop over triangles
    for f = 1:size(F, 1)
        if dblA(f) < min_dblarea
            % Find shortest edge
            E = [F(f,iedge(1,:));F(f,iedge(2,:));F(f,iedge(3,:))];
            Vdif = V(E(:,1),:)-V(E(:,2),:);
            d = [norm(Vdif(1,:)) norm(Vdif(2,:)) norm(Vdif(3,:))];
            [dmin,imin] = min(d);
            % Collapse min edge: i-->j
            i = iedge(imin,1);
            j = iedge(imin,2);
            FIM(F(f, i)) = FIM(F(f, j));
            num_edge_collapses = num_edge_collapses + 1;
        end
    end
    
    % Reindex faces
    rF = F;
    
    % Loop over triangles
    for f = 1:size(rF, 1)
        for i = 1:size(rF, 2)
            rF(f, i) = FIM(rF(f, i));
        end
    end
    
    indexv = (FIM)';
    FF = [];
    num_face_collapses = 0;
    
    % Only keep uncollapsed faces
    for f = 1:size(rF, 1)
        collapsed = false;
        
        % Check if any indices are the same
        for i = 1:size(rF, 2)
            for j = i+1:size(rF, 2)
                if rF(f, i) == rF(f, j)
                    collapsed = true;
                    num_face_collapses = num_face_collapses + 1;
                    break;
                end
            end
            if collapsed
                break;
            end
        end
        
        if ~collapsed
            FF = [FF; rF(f, :)];
        end
    end
    
    if num_edge_collapses == 0
        % There must have been a "collapsed edge" in the input
        assert(num_face_collapses == 0);
        % Base case
        return;
    end
    

  
end
