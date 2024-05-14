  
function [midpoint] = smooth_midpoint(point_1,point_2,normal_1,normal_2)
   
    v1 = point_1;
    v2 = point_2;
    n1 = normal_1;
    n2 = normal_2;
    
    %% Normals in the edges of the base are corrected
    nedges = size(point_1);
    for i = 1:nedges
        if (point_1(i,3) == 0) && (point_2(i,3) == 0)
            normal_1(i,3) = 0; normal_1(i,:) = normal_1(i,:)/norm(normal_1(i,:));
            normal_2(i,3) = 0; normal_2(i,:) = normal_2(i,:)/norm(normal_2(i,:));
        end
    end
    
    %% Bezier control points 
    %according to Lamberts, "Interpolation of curves and surfaces"
    direction_vector = point_2 - point_1;
    normal_midpoint = (normal_1 + normal_2)/2;
    normal_plane = cross(normal_midpoint,direction_vector,2);
    tangent_1 = cross(normal_plane,normal_1,2);
    tangent_2 = cross(normal_plane,normal_2,2);
    
    %Normalize the tangent vectors
    tangent_1 = tangent_1 ./ vecnorm(tangent_1,2,2);
    tangent_2 = tangent_2 ./ vecnorm(tangent_2,2,2);
    
    %Angles between the tangent and direction vector
    norm_direction_vector = vecnorm(direction_vector,2,2);
    theta = abs(180/pi*acos(dot(direction_vector, tangent_1,2))./norm_direction_vector);
    phi = abs(180/pi*acos(dot(direction_vector, tangent_2,2))./norm_direction_vector);
    
    %Computing alpha and beta:
    if (abs(theta) < 60)
        alpha = norm_direction_vector.^2./(3*dot(direction_vector,tangent_1,2));
    else
        alpha = norm_direction_vector*(2/3);
    end
    if (abs(phi) < 60)
        beta = norm_direction_vector.^2./(3*dot(direction_vector,tangent_2,2));
    else
        beta = norm_direction_vector*(2/3);
    end
    
    % Calculate the control point 1
    control_point_1 = point_1 + alpha .* tangent_1;
    
    % Calculate the control point 2
    control_point_2 = point_2 - beta .* tangent_2;

    %% Bezier curves.
    % https://www.youtube.com/watch?v=8Jb2f9R6br8&ab_channel=Prot04m
    numControl = 4;
    nPoints = 0.01;

    %P1, P4 are endpoints of the curve
    %P2, P3 are control points
    P1 = point_1;
    P2 = control_point_1;
    P3 = control_point_2;
    P4 = point_2;

    x = [P1(:,1) P2(:,1) P3(:,1) P4(:,1)];
    y = [P1(:,2) P2(:,2) P3(:,2) P4(:,2)];
    z = [P1(:,3) P2(:,3) P3(:,3) P4(:,3)];

    n =  numControl - 1;
    i = 0:n;
    coeff = factorial(n)./(factorial(i).*factorial(n-i));
    t = 0:nPoints:1;

    for j = 1:numControl
        b(j,:) = coeff(j)*t.^i(j).*(1 - t).^(n-i(j));
    end

    % Number of edges I will subdivide
    x_bezier = zeros(nedges(1),numel(t));
    y_bezier = zeros(nedges(1),numel(t));
    z_bezier = zeros(nedges(1),numel(t));
    for j = 1:numControl
        for k = 1:nedges(1)
            x_bezier(k,:) = b(j,:)*x(k,j) + x_bezier(k,:);
            y_bezier(k,:) = b(j,:)*y(k,j) + y_bezier(k,:);
            z_bezier(k,:) = b(j,:)*z(k,j) + z_bezier(k,:);
        end
    end

    %% Intersection between average normal at midpoint, and Bezier curve
    %%%%%% Run over the each of the edges
    tv = [transpose(t),transpose(t),transpose(t)];
    midpoint = [];
    normal_midpoint_curve2 = [];
    x_normal = zeros(nedges(1),numel(t));
    y_normal = zeros(nedges(1),numel(t));
    z_normal = zeros(nedges(1),numel(t));
    for k = 1:nedges(1)
        normal_midpoint_curve =(P1(k,:) + direction_vector(k,:)/2) + tv.*normal_midpoint(k,:);
        bezier_curve = [transpose(x_bezier(k,:)),transpose(y_bezier(k,:)),transpose(z_bezier(k,:))];
        %Matrix of distances from all points of the curves
        param_len = length(t);
        distances = zeros(param_len,param_len);
        for i = 1:param_len
            for j = 1:param_len
                distances(i,j) = sum((normal_midpoint_curve(i,:)-bezier_curve(j,:)).^2);
            end
        end
        distances_min = min(min(distances));
        [imid,jmid] = find(distances_min == distances);
        midpoint = [midpoint; [x_bezier(k,jmid),y_bezier(k,jmid),z_bezier(k,jmid)]];
        x_normal(k,:) = normal_midpoint_curve(:,1);
        y_normal(k,:) = normal_midpoint_curve(:,2);
        z_normal(k,:) = normal_midpoint_curve(:,3);
    end
    
%     tangent_1r = point_1 + tangent_1;
%     tangent_2r = point_2 + tangent_2;
%     figure(3); plot3(x,y,z,'o','MarkerFaceColor', 'k'); hold on
%     xlim([-1.00 1.50])
%     ylim([-1.000 1.500])
%     zlim([-1.00 1.50])
%     view([-89.80 0.36])
%     for k = 1:nedges(1)
%         plot3(midpoint(k,1),midpoint(k,2),midpoint(k,3),'o','MarkerFaceColor', 'g'); hold on
%         plot3(x_bezier(k,:),y_bezier(k,:),z_bezier(k,:), 'r'); hold on
%         plot3(x_normal(k,:),y_normal(k,:),z_normal(k,:), 'r'); hold on
%         plot3([point_1(k,1) tangent_1r(k,1)],[point_1(k,2) tangent_1r(k,2)],[point_1(k,3) tangent_1r(k,3)],'b-^', 'LineWidth',1); hold on
%         plot3([point_2(k,1) tangent_2r(k,1)],[point_2(k,2) tangent_2r(k,2)],[point_2(k,3) tangent_2r(k,3)],'b-^', 'LineWidth',1); hold on
%     end
%     hold off
end

   
  