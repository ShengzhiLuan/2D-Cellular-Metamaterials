%  Measurement of Features for 2D Cellular Metamaterials
%  Johns Hopkins University
%  Shengzhi Luan
%  02.08.2023
% =========================================================================
function  [ Feature ] = Feature_Measurement ...
                            ( Seed , Triangle , Center , p_Geometry , ...
                                t_Geometry , Box_Crop , Material_Property )
%  ------------------------------------------------------------------------
%  Input:  Seed -- Seed matrix ( 1st row is x-coordinate,
%                                2nd row is y-coordinate,
%                                3rd row is weight. ).
%          Triangle -- Laguerre triangle by anticlockwise.
%          Center -- Laguerre center of each Laguerre triangle:
%                                    ( 1st row is x-coordinate, 
%                                      2nd row is y-coordinate. ).
%          p_Geometry -- Node matrix ( 1st row is x-coordinate, 
%                                      2nd row is y-coordinate. ).
%          t_Geometry -- Connection matrix ( 1st row is start node ID,
%                                            2nd row is end node ID. ).
%          Box_Crop -- Cropping box ( [ x_min , x_max , y_min , y_max ] ).
%          Material_Property -- Material ( [Thickness , Stiffness ]).
%  Output: Feature -- Features of lattice material:
%                     1:        Relative density.
%                     2:        Number of polygon cells.
%                     3-6:      Cell edge number.
%                     7-10:     Cell area.
%                     11-14:    Cell compactness.
%                     15-18:    Cell eccentricity.
%                     19-22:    Strut length.
%                     23-26:    Strut orientation.
%                     27-30:    Nodal connectivity.
%                     31-34:    Neighboring seed distance.
%                     35-38:    Neighboring seed angle sine.
%                     39-42:    Neighboring seed angle cosine.
%  ------------------------------------------------------------------------
[ Density , Ligament_Length , Ligament_Angle , ...
            Node_Connection_Number ] = Feature_Lattice ...
                ( p_Geometry , t_Geometry , Material_Property , Box_Crop );
[ Seed_Near_Distance , Sin_Angle_Seed , Cos_Angle_Seed ] = ...
                    Feature_Center ( Seed , Triangle , Center , Box_Crop );
[ Cell_Number , Edge_Number , Area , Compactness , Eccentricity ] = ...
                      Feature_Cell ( Seed , Triangle , Center , Box_Crop );
%  ------------------------------------------------------------------------
Feature(1) = Density;
Feature(2) = Cell_Number;
Feature(3:6) = [mean(Edge_Number), ...
                var(Edge_Number), ...
                skewness(Edge_Number), ...
                kurtosis(Edge_Number)];
Feature(7:10) = [mean(Area), ...
                 var(Area), ...
                 skewness(Area), ...
                 kurtosis(Area)];
Feature(11:14) = [mean(Compactness), ...
                  var(Compactness), ...
                  skewness(Compactness), ...
                  kurtosis(Compactness)];
Feature(15:18) = [mean(Eccentricity), ...
                  var(Eccentricity), ...
                  skewness(Eccentricity), ...
                  kurtosis(Eccentricity)];
Feature(19:22) = [mean(Ligament_Length), ...
                  var(Ligament_Length), ...
                  skewness(Ligament_Length), ...
                  kurtosis(Ligament_Length)];
Feature(23:26) = [mean(Ligament_Angle), ...
                  var(Ligament_Angle), ...
                  skewness(Ligament_Angle), ...
                  kurtosis(Ligament_Angle)];
Feature(27:30) = [mean(Node_Connection_Number), ...
                  var(Node_Connection_Number), ...
                  skewness(Node_Connection_Number), ...
                  kurtosis(Node_Connection_Number)];
Feature(31:34) = [mean(Seed_Near_Distance), ...
                  var(Seed_Near_Distance), ...
                  skewness(Seed_Near_Distance), ...
                  kurtosis(Seed_Near_Distance)];
Feature(35:38) = [mean(Sin_Angle_Seed), ...
                  var(Sin_Angle_Seed), ...
                  skewness(Sin_Angle_Seed), ...
                  kurtosis(Sin_Angle_Seed)];
Feature(39:42) = [mean(Cos_Angle_Seed), ...
                  var(Cos_Angle_Seed), ...
                  skewness(Cos_Angle_Seed), ...
                  kurtosis(Cos_Angle_Seed)];
%  ------------------------------------------------------------------------
end
% =========================================================================




%  ============================== Subroutine ==============================
function  [ Density , Ligament_Length , Ligament_Angle , ...
            Node_Connection_Number ] = Feature_Lattice ...
                 ( p_Geometry , t_Geometry , Material_Property , Box_Crop )
%  ------------------------------------------------------------------------
%  Input:  p_Geometry -- Node matrix ( 1st row is x-coordinate, 
%                                      2nd row is y-coordinate. ).
%          t_Geometry -- Connection matrix ( 1st row is start node ID,
%                                            2nd row is end node ID. ).
%          Material_Property -- Vector of material property 
%                               [ Thickness of ligament ].
%          Box_Crop -- Cropping box ( [ x_min , x_max , y_min , y_max ] ).
%  Output: Density -- Relative density.
%          Ligament_Length -- Length of each ligament.
%          Ligament_Angle -- Angle of each ligament.
%          Node_Connection_Number -- Number of connection at each node.
%  ------------------------------------------------------------------------
p = p_Geometry;
t = t_Geometry;
Lx = Box_Crop(2)-Box_Crop(1);
Ly = Box_Crop(4)-Box_Crop(3);
Thickness = Material_Property(1);
%  ------------------------------------------------------------------------
Node_Flag = zeros(1,size(p,2));
for i = 1:1:size(p,2)
    if p(1,i) >= Box_Crop(1) && p(1,i) <= Box_Crop(2)
        if p(2,i) >= Box_Crop(3) && p(2,i) <= Box_Crop(4)
            Node_Flag(i) = 1;
        end
    end
end
%  ------------------------------------------------------------------------
count = 0;
for i = 1:1:size(p,2)
    if Node_Flag(i) == 1
        count = count + 1;
        Node_Connection_Number(count) = length(find(t==i));
    end
end
%  ------------------------------------------------------------------------
count = 0;
for i = 1:1:size(t,2)
    if sum(Node_Flag(t(:,i)')) ~= 0
        count = count + 1;
        Ligament_Length(count) = norm(p(:,t(1,i))-p(:,t(2,i)));
    end
end
%  ------------------------------------------------------------------------
[ p_Crop , t_Crop ] = Box_Cropping ( p_Geometry , t_Geometry , Box_Crop );
p = p_Crop;
t = t_Crop;
for i = 1:1:size(t,2)
    Ligament_Vector = p(:,t(1,i))-p(:,t(2,i));
    Ligament_Crop_Length(i) = norm(Ligament_Vector);
    Angle_1 = acos([0,1]*Ligament_Vector/norm(Ligament_Vector))/pi*180;
    Angle_2 = acos([0,1]*-Ligament_Vector/norm(-Ligament_Vector))/pi*180;
    Ligament_Angle(i) = min([Angle_1,Angle_2]);
end
Density = sum(Ligament_Crop_Length)*Thickness/Lx/Ly;
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Seed_Near_Distance , Sin_Angle_Seed , ...
            Cos_Angle_Seed ] = Feature_Center ...
                                    ( Seed , Triangle , Center , Box_Crop )
%  ------------------------------------------------------------------------
%  Input:  Seed -- Seed matrix ( 1st row is x-coordinate, 
%                                2nd row is y-coordinate,
%                                3rd row is weight. ).
%          Triangle -- Laguerre triangle by anticlockwise.
%          Center -- Laguerre center of each Laguerre triangle:
%                                    ( 1st row is x-coordinate, 
%                                      2nd row is y-coordinate. ).
%          Box_Crop -- Cropping box ( [ x_min , x_max , y_min , y_max ] ).
%  Output: Seed_Near_Distance -- Nearest distance around neighbor seeds.
%          Cos_Angle_Seed -- Average cosine angle around neighbor seeds.
%          Sin_Angle_Seed -- Average sine angle around neighbor seeds.
%  ------------------------------------------------------------------------
p = Seed([1,2],:);
for i = 1:1:size(Seed,2)
    if p(1,i) >= Box_Crop(1) && p(1,i) <= Box_Crop(2)
        if p(2,i) >= Box_Crop(3) && p(2,i) <= Box_Crop(4)
            [~,ID] = find(Triangle==i);
            Cell_Matrix = uniquetol(Center(:,ID)',1E-3,'Byrows',true)';
            Cell_Center = [mean(Cell_Matrix(1,:))
                           mean(Cell_Matrix(2,:))];
            p(:,i) = Cell_Center;
        end
    end
end
%  ------------------------------------------------------------------------
for i = 1:1:size(Triangle,2)
    t_Original(:,3*(i-1)+1) = [Triangle(1,i);Triangle(2,i)];
    t_Original(:,3*(i-1)+2) = [Triangle(2,i);Triangle(3,i)];
    t_Original(:,3*(i-1)+3) = [Triangle(3,i);Triangle(1,i)];
end
for i = 1:1:size(t_Original,2)
    if t_Original(1,i) < t_Original(2,i)
        t_Duplicated(1,i) = t_Original(1,i);
        t_Duplicated(2,i) = t_Original(2,i);
    end
    if t_Original(1,i) > t_Original(2,i)
        t_Duplicated(1,i) = t_Original(2,i);
        t_Duplicated(2,i) = t_Original(1,i);
    end
end
t = unique(t_Duplicated','rows')';
%  ------------------------------------------------------------------------
count = 0;
for i = 1:1:size(p,2)
    if p(1,i) >= Box_Crop(1) && p(1,i) <= Box_Crop(2)
        if p(2,i) >= Box_Crop(3) && p(2,i) <= Box_Crop(4)
            [I,J] = find(t==i);
            if length(J) ~= 0
                count = count + 1;
                for j = 1:1:length(J)
                    Dy = p(2,t(3-I(j),J(j)))-p(2,t(I(j),J(j)));
                    Dx = p(1,t(3-I(j),J(j)))-p(1,t(I(j),J(j)));
                    Seed_Distance(j) = sqrt(Dx^2+Dy^2);
                    Seed_Angle(j) = atan2(Dy,Dx);
                end
                Seed_Near_Distance(count) = min(Seed_Distance);
                Sin_Angle_Seed(count) = mean(sin(Seed_Angle));
                Cos_Angle_Seed(count) = mean(cos(Seed_Angle));
                clear Seed_Distance Seed_Angle
            end
        end
    end
end
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Cell_Number , Edge_Number , Area , Compactness , ...
            Eccentricity ] = Feature_Cell ...
                                    ( Seed , Triangle , Center , Box_Crop )
%  ------------------------------------------------------------------------
%  Input:  Seed -- Seed matrix ( 1st row is x-coordinate, 
%                                2nd row is y-coordinate,
%                                3rd row is weight. ).
%          Triangle -- Laguerre triangle by anticlockwise.
%          Center -- Laguerre center of each Laguerre triangle:
%                                    ( 1st row is x-coordinate, 
%                                      2nd row is y-coordinate. ).
%          Box_Crop -- Cropping box ( [ x_min , x_max , y_min , y_max ] ).
%  Output: Cell_Number -- Number of cell.
%          Edge_Number -- Edge number of each cell.
%          Area -- Area of each cell.
%          Compactness -- Compactness of each cell.
%          Roundness -- Roundness of each cell.
%          Eccentricity -- Eccentricity of each cell.
%  ------------------------------------------------------------------------
p = Seed(1:2,:);
count = 0;
for i = 1:1:size(p,2)
    if p(1,i) >= Box_Crop(1) && p(1,i) <= Box_Crop(2)
        if p(2,i) >= Box_Crop(3) && p(2,i) <= Box_Crop(4)
            [~,ID] = find(Triangle==i);
            if length(ID) ~= 0
                Cell_Matrix = uniquetol(Center(:,ID)',1E-3,'Byrows',true)';
                Cell_Center = [mean(Cell_Matrix(1,:))
                               mean(Cell_Matrix(2,:))];
                Node_Angle = atan2(Cell_Matrix(2,:)-Cell_Center(2),...
                                   Cell_Matrix(1,:)-Cell_Center(1));
                [~,Cell_Order] = sort(Node_Angle);
                Cell_Matrix = Cell_Matrix(:,Cell_Order);
                %  --------------------------------------------------------
                if size(Cell_Matrix,2) >= 3
                    count = count + 1;
                    Area(count) = 0;
                    Perimeter = 0;
                    for j = 1:1:size(Cell_Matrix,2)
                        Cell = [Cell_Matrix,Cell_Matrix(:,1)];
                        Cell_Area = 1/2*det([Cell(1,j),Cell(2,j),1
                                         Cell(1,j+1),Cell(2,j+1),1
                                   Cell_Center(1),Cell_Center(2),1]);
                        Area(count) = Area(count)+Cell_Area;
                        Perimeter = Perimeter+norm(Cell(:,j)-Cell(:,j+1));
                    end
                    [ Axis_Length ] = Major_Minor_Axis ( Cell_Matrix );
                    Compactness(count) = 4*pi*Area(count)/Perimeter^2;
                    Eccentricity(count) = Axis_Length(1)/Axis_Length(2);
                    Edge_Number(count) = size(Cell_Matrix,2);
                end
            end
        end
    end
end
%  ------------------------------------------------------------------------
Cell_Number = count;
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Axis_Length ] = Major_Minor_Axis ( p_Input )
%  ------------------------------------------------------------------------
%  Input:  p_Input -- Node matrix ( 1st row is x-coordinate;
%                                   2nd row is y-coordinate. ).
%  Output: Aixs_Length -- Length of minor and major axis.
%  ------------------------------------------------------------------------
p = p_Input;
%  ------------------------------------------------------------------------
Hull = convhull(p(1,:),p(2,:));
CH = p(:,Hull);
E = diff(CH,1,2);
T = atan2(E(2,:),E(1,:));
T = unique(mod(T,pi/2));
%  ------------------------------------------------------------------------
R = cos( reshape(repmat(T,2,2),2*length(T),2) ...
       + repmat([0 -pi ; pi 0]/2,length(T),1));
RCH = R*CH;
bsize = max(RCH,[],2) - min(RCH,[],2);
area  = prod(reshape(bsize,2,length(bsize)/2));
[~,i] = min(area);
%  ------------------------------------------------------------------------
Rf    = R(2*i+[-1 0],:);
bound = Rf * CH;
bmin  = min(bound,[],2);
bmax  = max(bound,[],2);
Rf = Rf';
bb(:,4) = bmax(1)*Rf(:,1) + bmin(2)*Rf(:,2);
bb(:,1) = bmin(1)*Rf(:,1) + bmin(2)*Rf(:,2);
bb(:,2) = bmin(1)*Rf(:,1) + bmax(2)*Rf(:,2);
bb(:,3) = bmax(1)*Rf(:,1) + bmax(2)*Rf(:,2);
%  ------------------------------------------------------------------------
bb = [bb,bb(:,1)];
for i = 1:1:size(bb,2)-1
    Length(i) = norm(bb(:,i)-bb(:,i+1));
end
Axis_Length = [min(Length),max(Length)];
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ p_Output , t_Output ] = Box_Cropping ...
                                    ( p_Input , t_Input , Box_Crop )
%  ------------------------------------------------------------------------
%  Input:  p_Input -- Nodal matrix ( 1st row is x-coordinate,
%                                    2nd row is y-coordinate,
%                                    3rd row is z-coordinate. ).
%          t_Input -- Connection matrix ( 1st row is start node ID,
%                                         2nd row is end node ID. ).
%          Box -- Cropping box ([x_min,x_max,y_min,y_max]).
%  Output: p_Output -- Cropped node matrix with same format.
%          t_Output -- Cropped connect matrix with same format.
%  ------------------------------------------------------------------------
Ligament_Ann_Coe = 0.15;
p = p_Input;
t = t_Input;
for i = 1:1:size(t,2)
    Length(i) = norm(p(:,t(1,i))-p(:,t(2,i)));
end
Length_Mean = sum(Length)/length(Length);
Tol = 1E-3*Length_Mean;
x1 = Box_Crop(1);
x2 = Box_Crop(2);
y1 = Box_Crop(3);
y2 = Box_Crop(4);
[ a , b , c ] = Line_Equation ( [x1,y1] , [x2,y1] );
Box_Line_1 = [a,b,c];
[ a , b , c ] = Line_Equation ( [x2,y1] , [x2,y2] );
Box_Line_2 = [a,b,c];
[ a , b , c ] = Line_Equation ( [x2,y2] , [x1,y2] );
Box_Line_3 = [a,b,c];
[ a , b , c ] = Line_Equation ( [x1,y2] , [x1,y1] );
Box_Line_4 = [a,b,c];
%  ------------------------------------------------------------------------
Node_Flag = zeros(1,size(p,2));
x_Box = [x1,x2,x2,x1,x1];
y_Box = [y1,y1,y2,y2,y1];
[In_Index,On_Index] = inpolygon(p(1,:),p(2,:),x_Box,y_Box);
for i = 1:1:size(p,2)
    if In_Index(i) == 0
        Node_Flag(i) = 0;
    elseif On_Index(i) == 1
        Node_Flag(i) = 1;
    else
        Node_Flag(i) = 2;
    end
end
%  ------------------------------------------------------------------------
Element_Flag = zeros(1,size(t,2));
count = 0;
for i = 1:1:size(t,2)
    Flag_Vector = sort([Node_Flag(t(1,i)),Node_Flag(t(2,i))]);
    if Flag_Vector(1) == 0 && Flag_Vector(2) == 2
        Element_Flag(i) = 1;
        ID = find(Node_Flag(t(:,i))==0);
        [ a , b , c ] = Line_Equation ( p(:,t(1,i))' , p(:,t(2,i))' );
        [ x_1 , y_1 ] = Intersect_Point ( [a,b,c] , Box_Line_1 );
        [ x_2 , y_2 ] = Intersect_Point ( [a,b,c] , Box_Line_2 );
        [ x_3 , y_3 ] = Intersect_Point ( [a,b,c] , Box_Line_3 );
        [ x_4 , y_4 ] = Intersect_Point ( [a,b,c] , Box_Line_4 );
        Inter_Matrix = [x_1,x_2,x_3,x_4;y_1,y_2,y_3,y_4];
        for j = 1:1:4
            x_Inter = Inter_Matrix(1,j);
            y_Inter = Inter_Matrix(2,j);
            x_Sign = sign((p(1,t(1,i))-x_Inter)*(p(1,t(2,i))-x_Inter));
            y_Sign = sign((p(2,t(1,i))-y_Inter)*(p(2,t(2,i))-y_Inter));
            Inter_Flag(:,j) = [j;x_Sign+y_Sign
                               norm(p(:,t(3-ID,i))-Inter_Matrix(:,j))];
        end
        Inter_Flag_Mid = Inter_Flag(:,find(Inter_Flag(2,:)<2));
        [~,Inter_ID] = min(Inter_Flag_Mid(3,:));
        Inter_ID = Inter_Flag_Mid(1,Inter_ID);
        count = count + 1;
        t(ID,i) = count + size(p,2);
        p_New(:,count) = Inter_Matrix(:,Inter_ID);
        if Inter_ID == 1
            p_New(2,count) = y1;
        end
        if Inter_ID == 2
            p_New(1,count) = x2;
        end
        if Inter_ID == 3
            p_New(2,count) = y2;
        end
        if Inter_ID == 4
            p_New(1,count) = x1;
        end
    end
end
%  ------------------------------------------------------------------------
p = [p,p_New];
Delete_Node = find(Node_Flag==0);
[ t_New ] = Node_Deletion ( p , t , Delete_Node );
p(:,Delete_Node) = [];
%  ------------------------------------------------------------------------
p_Input = p;
t_Input = t_New;
[ p_Output , t_Output ] = Ligament_Annealing ...
                                  ( p_Input , t_Input , Ligament_Ann_Coe );
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ x , y ] = Intersect_Point ( Line_1 , Line_2 )
%  ------------------------------------------------------------------------
%  Input:  Line_1 -- Coefficients of Line 1 (a,b,c).
%          Line_1 -- Coefficients of Line 2 (a,b,c).
%  Output: [x,y] -- Coordinates of intesection point.
%  ------------------------------------------------------------------------
a1 = Line_1(1);
b1 = Line_1(2);
c1 = Line_1(3);
a2 = Line_2(1);
b2 = Line_2(2);
c2 = Line_2(3);
%  ------------------------------------------------------------------------
if a1*b2 ~= a2*b1
    x = (b1*c2-b2*c1)/(a1*b2-a2*b1);
    y = (a2*c1-a1*c2)/(a1*b2-a2*b1);
end
if a1*b2 == a2*b1
    x = Inf;
    y = Inf;
end
%  ------------------------------------------------------------------------
end
% =========================================================================




%  ============================== Subroutine ==============================
function  [ a , b , c ] = Line_Equation ( Point_A , Point_B )
%  ------------------------------------------------------------------------
%  Input:  Point_A -- Coordinates of point A (x,y).
%          Point_B -- Coordinates of point B (x,y).
%  Output: [a,b,c] -- Coefficients in equation: a*x+b*y+c=0.
%  ------------------------------------------------------------------------
x1 = Point_A(1);
y1 = Point_A(2);
x2 = Point_B(1);
y2 = Point_B(2);
%  ------------------------------------------------------------------------
if x1 ~= x2
    k = (y2-y1)/(x2-x1);
    l = y1-k*x1;
    a = k;
    b = -1;
    c = l;
end
if x1 == x2
    a = 1;
    b = 0;
    c = -x1;
end
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ M_Output ] = Node_Deletion ( p_Input , M_Input , Node_ID )
%  ------------------------------------------------------------------------
%  Input:  p_Input -- Node matrix.
%          M_Input -- Original input matrix (pair in column).
%          Node_ID -- Node Deletion ID.
%  Output: M_Output -- Modified output matrix (pair in column).
%  ------------------------------------------------------------------------
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