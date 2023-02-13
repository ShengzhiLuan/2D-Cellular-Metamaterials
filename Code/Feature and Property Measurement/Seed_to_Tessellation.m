%  Generation of Tessellation by Seed & Weight
%  Johns Hopkins University
%  Shengzhi Luan
%  02.08.2023
% =========================================================================
function  [ Seed_Output , p_Output , t_Output , Triangle , Center ] = ...
                                        Seed_to_Tessellation ( Seed_Input )
%  ------------------------------------------------------------------------
%  Input:  Seed_Input -- Seed matrix ( 1st row is x-coordinate,
%                                      2nd row is y-coordinate,
%                                      3rd row is weight. ).
%  Output: Seed_Output -- Seed matrix ( 1st row is x-coordinate,
%                                       2nd row is y-coordinate,
%                                       3rd row is weight. ).
%          p_Output -- Node matrix ( 1st row is x-coordinate,
%                                    2nd row is y-coordinate. ).
%          t_Output -- Connection matrix ( 1st row is start node ID,
%                                          2nd row is end node ID. ).
%          Triangle -- Laguerre triangle in anticlockwise direction.
%          Center -- Laguerre center ( 1st row is x-coordinate,
%                                      2nd row is y-coordinate. ).
%  ------------------------------------------------------------------------
Ligament_Anneal_Coefficient = 0.15;
Seed_Anneal_Coefficient = 0.15;
Seed = Seed_Input;
[ Seed , Triangle , Center ] = Laguerre_Triangulation ...
                                        ( Seed , Seed_Anneal_Coefficient );
[ p_Frame , t_Frame ] = Laguerre_Tessellation ( Seed , Triangle , Center );
[ p_Frame , t_Frame ] = Ligament_Annealing ...
                       ( p_Frame , t_Frame , Ligament_Anneal_Coefficient );
%  ------------------------------------------------------------------------
for i = 1:1:size(Center,2)
    Center_to_Node = Center(:,i)-p_Frame;
    Distance = [];
    for j = 1:1:size(Center_to_Node,2)
        Distance(j) = norm(Center_to_Node(:,j));
    end
    [~,ID] = min(Distance);
    Center(:,i) = p_Frame(:,ID);
end
%  ------------------------------------------------------------------------
Seed_Output = Seed;
p_Output = p_Frame;
t_Output = t_Frame;
%  ------------------------------------------------------------------------
end
% =========================================================================




%  ============================== Subroutine ==============================
function  [ Seed_Output , Triangle , Center ] = Laguerre_Triangulation ...
                                   ( Seed_Input , Seed_Anneal_Coefficient )
%  ------------------------------------------------------------------------
%  Input:  Seed_Input -- Seed input matrix ( 1st row is x-coordinate,
%                                            2nd row is y-coordinate,
%                                            3rd row is weight. ).
%          Seed_Anneal_Coefficient -- Seed annealing coefficient.
%  Output: Seed_Output -- Seed output matrix ( 1st row is x-coordinate,
%                                              2nd row is y-coordinate,
%                                              3rd row is weight. ).
%          Triangle -- Laguerre triangle in anticlockwise direction.
%          Center -- Laguerre center ( 1st row is x-coordinate,
%                                      2nd row is y-coordinate. ).
%  ------------------------------------------------------------------------
[ Triangle , Center ] = Power_Triangulation ( Seed_Input );
[ Seed_Output ] = Seed_Annealing ...
              ( Seed_Input , Triangle , Center , Seed_Anneal_Coefficient );
[ Triangle , Center ] = Power_Triangulation ( Seed_Output );
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ p_Output , t_Output ] = Laguerre_Tessellation ...
                                               ( Seed , Triangle , Center )
%  ------------------------------------------------------------------------
%  Input:  Seed -- Seed matrix ( 1st row is x-coordinate,
%                                2nd row is y-coordinate,
%                                3rd row is weight. ).
%          Triangle -- Laguerre triangle in anticlockwise direction.
%          Center -- Laguerre center ( 1st row is x-coordinate,
%                                      2nd row is y-coordinate. ).
%  Output: p_Output -- Node matrix ( 1st row is x-coordinate,
%                                    2nd row is y-coordinate. ).
%          t_Output -- Connection matrix ( 1st row is start node ID,
%                                          2nd row is end node ID. ).
%  ------------------------------------------------------------------------
Seed_Box = [min(Seed(1,:)),max(Seed(1,:)),min(Seed(2,:)),max(Seed(2,:))];
Dimension = [Seed_Box(1:2),Seed_Box(2)-Seed_Box(1), ...
             Seed_Box(3:4),Seed_Box(4)-Seed_Box(3)];
%  ------------------------------------------------------------------------
Tri = Triangle';
Cen = Center';
Tri_Vector = repmat((1:size(Tri,1))',1,3);
Tri_Matrix = sparse(Tri,Tri(:,[3 1 2]),Tri_Vector, ...
                    size(Seed,2),size(Seed,2));
Inter_Edge = (Tri_Matrix & Tri_Matrix').*Tri_Matrix; 
Outer_Edge = xor(Tri_Matrix, Tri_Matrix').*Tri_Matrix;
[~,~,v] = find(triu(Inter_Edge));
[~,~,vv] = find(triu(Inter_Edge'));
vector_x = [Cen(v,1) Cen(vv,1)]';
vector_y = [Cen(v,2) Cen(vv,2)]';
[i,j,z] = find(Outer_Edge);
dx = Seed(1,j)' - Seed(1,i)';
dy = Seed(2,j)' - Seed(2,i)';
rx = max(Seed(1,:)')-min(Seed(2,:)'); 
ry = max(Seed(1,:)')-min(Seed(2,:)');
cx = (max(Seed(1,:)')+min(Seed(1,:)'))/2 - Cen(z,1); 
cy = (max(Seed(2,:)')+min(Seed(2,:)'))/2 - Cen(z,2);
nm = sqrt(rx.*rx + ry.*ry) + sqrt(cx.*cx + cy.*cy);
scale = nm./sqrt((dx.*dx+dy.*dy));
edge_x = [Cen(z,1) Cen(z,1)-dy.*scale]';
edge_y = [Cen(z,2) Cen(z,2)+dx.*scale]';
Vector_X = [vector_x edge_x];
Vector_Y = [vector_y edge_y];
%  ------------------------------------------------------------------------
Vector = [Vector_X;Vector_Y];
count = 0;
Delete_ID = [];
for i = 1:1:size(Vector,2)
    if length(find(Vector(1:2,i)<Dimension(1)-0.01*Dimension(3)))~=0 || ...
       length(find(Vector(1:2,i)>Dimension(2)+0.01*Dimension(3)))~=0 || ...
       length(find(Vector(3:4,i)<Dimension(4)-0.01*Dimension(6)))~=0 || ...
       length(find(Vector(3:4,i)>Dimension(5)+0.01*Dimension(6)))~=0
        count = count + 1;
        Delete_ID(count) = i;
    end
end
Vector(:,Delete_ID) = [];
Vector_X = Vector(1:2,:);
Vector_Y = Vector(3:4,:);
%  ------------------------------------------------------------------------
for i = 1:1:size(Vector_X,2)
    p(:,2*i-1) = [Vector_X(1,i);Vector_Y(1,i)];
    p(:,2*i) = [Vector_X(2,i);Vector_Y(2,i)];
    t(:,i) = [2*i-1;2*i];
end
[p_Unique,~,Index] = uniquetol(p',1E-3,'Byrows',true);
p_Unique = p_Unique';
for i = 1:1:size(t,1)
    for j = 1:1:size(t,2)
        t_Updated(i,j) = Index(t(i,j));
    end
end
count = 0;
for i = 1:1:size(t_Updated,2)
    if t_Updated(1,i) ~= t_Updated(2,i)
        count = count + 1;
        t_Unique(:,count) = t_Updated(:,i);
    end
end
%  ------------------------------------------------------------------------
p_Output = p_Unique;
t_Output = t_Unique;
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ p_Output , t_Output ] = Ligament_Annealing ...
                        ( p_Input , t_Input , Ligament_Anneal_Coefficient )
%  ------------------------------------------------------------------------
%  Input:  p_Input -- Node matrix ( 1st row is x-coordinate,
%                                   2nd row is y-coordinate. ).
%          t_Input -- Connection matrix ( 1st row is start node ID,
%                                         2nd row is end node ID. ).
%          Ligament_Anneal_Coefficient -- Ligament annealing coefficient.
%  Output: p_Output -- Node matrix ( 1st row is x-coordinate,
%                                    2nd row is y-coordinate. ).
%          t_Output -- Connection matrix ( 1st row is start node ID,
%                                          2nd row is end node ID. ).
%  ------------------------------------------------------------------------
Annealing_Coefficient = Ligament_Anneal_Coefficient;
p = p_Input;
t = t_Input;
Number_of_Node = size(p,2);
p_Annealing = p;
t_Annealing = t;
%  ------------------------------------------------------------------------
for i = 1:1:size(t,2)
    Length(i) = norm(p(:,t(1,i))-p(:,t(2,i)));
end
Annealing_Element = find(Length <= Annealing_Coefficient*mean(Length));
%  ------------------------------------------------------------------------
if length(Annealing_Element) ~= 0
    for i = 1:1:length(Annealing_Element)
        Node_ID(2*i-1) = t(1,Annealing_Element(i));
        Node_ID(2*i) = t(2,Annealing_Element(i));
        p(:,Number_of_Node+i) = 1/2*(p(:,t(1,Annealing_Element(i))) + ...
                                     p(:,t(2,Annealing_Element(i))));
        t(t==t(1,Annealing_Element(i))) = Number_of_Node+i;
        t(t==t(2,Annealing_Element(i))) = Number_of_Node+i;
    end
    %  --------------------------------------------------------------------
    t(:,Annealing_Element) = [];
    Node_ID = unique(sort(Node_ID));
    [ t_Delete ] = Node_Deletion ( p , t , Node_ID );
    p(:,Node_ID) = [];
    %  --------------------------------------------------------------------
    p_Annealing = p;
    t_Annealing = t_Delete;
end
%  ------------------------------------------------------------------------
p_Output = p_Annealing;
t_Output = t_Annealing;
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Triangle , Center ] = Power_Triangulation ( Seed )
%  ------------------------------------------------------------------------
%  Input:  Seed -- Seed matrix ( 1st row is x-coordinate,
%                                2nd row is y-coordinate,
%                                3rd row is weight. ).
%  Output: Triangle -- Laguerre triangle in anticlockwise direction.
%          Center -- Laguerre center ( 1st row is x-coordinate,
%                                      2nd row is y-coordinate. ).
%  ------------------------------------------------------------------------
for i = 1:1:size(Seed,2)
    Seed_Lifted(1:2,i) = Seed(1:2,i);
    Seed_Lifted(3,i) = Seed(1,i)^2+Seed(2,i)^2-Seed(3,i)^2;
end
Cell = convhulln(Seed_Lifted')';
%  ------------------------------------------------------------------------
Seed_Lifted = Seed_Lifted';
Cell = Cell';
%  ------------------------------------------------------------------------
Model_Center = mean(Seed_Lifted,1);
for i = 1:1:size(Cell,1)
    vector = null(Seed_Lifted(Cell(i,1),:)-Seed_Lifted(Cell(i,2:3),:))';
    if size(vector,1) > 1
        Vector(i,:) = NaN;
    else
        Vector(i,:) = vector;
    end
    Triangle_Center(i,:) = mean(Seed_Lifted(Cell(i,:),:),1);
end
Dot_Flag = sum((Model_Center-Triangle_Center).*Vector,2);
Outer_ID = Dot_Flag < 0;
Vector(Outer_ID,:) = -1*Vector(Outer_ID,:);
Index = Vector(:,size(Cell,2)) > 0;
Index = Index';
%  ------------------------------------------------------------------------
Seed_Lifted = Seed_Lifted';
Cell = Cell';
Triangle = Cell(:,Index);
%  ------------------------------------------------------------------------
for i = 1:1:size(Triangle,2)
    A = Seed_Lifted(:,Triangle(1,i));
    B = Seed_Lifted(:,Triangle(2,i));
    C = Seed_Lifted(:,Triangle(3,i));
    N = cross(A,B)+cross(B,C)+cross(C,A);
    Center(i,:) = -0.5/N(3)*N(1:2);
end
Center = Center';
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Seed_Output ] = Seed_Annealing ...
            ( Seed_Input , Triangle , Center , Seed_Anneal_Coefficient )
%  ------------------------------------------------------------------------
%  Input:  Seed_Input -- Seed input matrix ( 1st row is x-coordinate,
%                                            2nd row is y-coordinate,
%                                            3rd row is weight. ).
%          Triangle -- Laguerre triangle in anticlockwise direction.
%          Center -- Laguerre center ( 1st row is x-coordinate,
%                                      2nd row is y-coordinate. ).
%          Seed_Anneal_Coefficient -- Seed annealing coefficient.
%  Output: Seed_Output -- Seed output matrix ( 1st row is x-coordinate,
%                                              2nd row is y-coordinate,
%                                              3rd row is weight. ).
%  ------------------------------------------------------------------------
Annealing_Coefficient = Seed_Anneal_Coefficient;
p = Seed_Input(1:2,:);
%  ------------------------------------------------------------------------
for i = 1:1:size(p,2)
    [~,ID] = find(Triangle==i);
    Cell_Matrix = uniquetol(Center(:,ID)',1E-3,'Byrows',true)';
    Cell_Center = [mean(Cell_Matrix(1,:))
                   mean(Cell_Matrix(2,:))];
    Node_Angle = atan2(Cell_Matrix(2,:)-Cell_Center(2),...
                       Cell_Matrix(1,:)-Cell_Center(1));
    [~,Cell_Order] = sort(Node_Angle);
    Cell_Matrix = Cell_Matrix(:,Cell_Order);
    Area(i) = 0;
    for j = 1:1:size(Cell_Matrix,2)
        Cell = [Cell_Matrix,Cell_Matrix(:,1)];
        Cell_Area = 1/2*det([Cell(1,j),Cell(2,j),1
                             Cell(1,j+1),Cell(2,j+1),1
                             Cell_Center(1),Cell_Center(2),1]);
        Area(i) = Area(i) + Cell_Area;
    end
end
%  ------------------------------------------------------------------------
Area_with_ID = [1:1:length(Area);Area];
Area_with_ID(:,isnan(Area_with_ID(2,:))) = [];
Box_Area = (range(p(1,:))*range(p(2,:)));
ID_Large = find(Area_with_ID(2,:)>Box_Area/size(p,2)*10);
Area_with_ID(:,ID_Large) = [];
ID_Small = Area_with_ID(1,find(Area_with_ID(2,:)< ...
                        mean(Area_with_ID(2,:))*Annealing_Coefficient));
%  ------------------------------------------------------------------------
Seed = Seed_Input;
Seed(:,ID_Small) = [];
Seed_Output = Seed;
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ M_Output ] = Node_Deletion ( p_Input , M_Input , Node_ID )
% -------------------------------------------------------------------------
%  Input:  p_Input -- Node matrix.
%          M_Input -- Original input matrix ( pair in column ).
%          Node_ID -- Node Deletion ID.
%  Output: M_Output -- Modified output matrix ( pair in column ).
% -------------------------------------------------------------------------
p = p_Input;
M = M_Input;
Node_ID_Mod = [0,Node_ID,size(p,2)+1];
%  ------------------------------------------------------------------------
for i = 1:1:size(M,1)
    for j = 1:1:size(M,2)
        for k = 1:1:length(Node_ID_Mod)-1
            if M(i,j) > Node_ID_Mod(k) && M(i,j) < Node_ID_Mod(k+1)
                M_New(i,j) = M(i,j) - (k-1);
            end
            if M(i,j) == Node_ID_Mod(k+1)
                M_New(i,j) = -1;
            end
        end
    end
end
%  ------------------------------------------------------------------------
[~,j] = find(M_New==-1);
M_New(:,j') = [];
M_Output = M_New;
%  ------------------------------------------------------------------------
end
%  ========================================================================