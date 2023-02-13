%  Measurement of Microstructural Feature and Macroscopic Property
%  Johns Hopkins University
%  Shengzhi Luan
%  02.08.2023
% =========================================================================
Sample_Number = 1646;
Density_Vector = [1,5,10,15,20]/100;
Box_Crop = [0.20,0.80,0.20,0.80];
V_Box = (Box_Crop(2)-Box_Crop(1))*(Box_Crop(4)-Box_Crop(3));
n_Element = 10;
Pause_Time = 5;
Data_Count = 0;
%  ------------------------------------------------------------------------
for i = 1:1:1
    fprintf('-- Measurement of %d / %d architected metamaterial ...\n', ...
            i,Sample_Number);
    cd('Seed_Data');
    File_ID = strcat('Seed_',num2str(i),'.txt');
    Seed = load(File_ID);
    cd ..
    [ Seed , p_Frame , t_Frame , ...
             Triangle , Center ] = Seed_to_Tessellation ( Seed );
    [ p_Crop , t_Crop ] = Box_Cropping ( p_Frame , t_Frame , Box_Crop );
    [ p_Lattice , t_Lattice ] = Beam_Discretization_Quadrature ...
                                ( p_Crop , t_Crop , n_Element );
    Sum_Length = 0;
    for j = 1:1:size(t_Crop,2)
        Sum_Length = Sum_Length + norm(p_Crop(:,t_Crop(1,j)) - ...
                                       p_Crop(:,t_Crop(2,j)));
    end
    for j = 1:1:length(Density_Vector)
        fprintf('  -- Analysis of %d / %d different density ...\n', ...
                j,length(Density_Vector));
        mkdir('Simulation_Folder');
        Data_Count = Data_Count + 1;
        V_Micro = V_Box*Density_Vector(j);
        Thickness = V_Micro/Sum_Length;
        V_Micro_Trial = 0;
        count = 0;
        while V_Micro_Trial < V_Micro && ...
              abs(V_Micro_Trial-V_Micro) >= 1E-5*V_Micro
            Polygon_Vector = [];
            count = count + 1;
            fprintf('Iteration %d: density accuracy as %f ...\n', ...
                    count, V_Micro_Trial/V_Micro);
            for ii = 1:1:size(t_Crop,2)
                Node_A = p_Crop(:,t_Crop(1,ii));
                Node_B = p_Crop(:,t_Crop(2,ii));
                [ Line ] = Polygon_Line ( Node_A , Node_B , Thickness );
                Polygon_Vector = [Polygon_Vector,Line];
            end
            [ Polygon_Lattice ] = Polygon_Union ( Polygon_Vector );
            V_Micro_Trial = Polygon_Lattice.area;
            Thickness = Thickness*V_Micro/V_Micro_Trial;
            if count >= 100
                break
            end
        end
        Material_Property = Thickness;
        [ Feature ] = Feature_Measurement ( Seed , Triangle , Center , ...
                        p_Frame , t_Frame , Box_Crop , Material_Property );
        Feature(1) = Density_Vector(j);
        copyfile Abaqus_Input.m Simulation_Folder;
        cd('Simulation_Folder');
        Bottom = find(p_Lattice(2,:)<=Box_Crop(3)+ ...
                      0.125*(Box_Crop(4)-Box_Crop(3)));
        Top = find(p_Lattice(2,:)>=Box_Crop(4)- ...
                   0.125*(Box_Crop(4)-Box_Crop(3)));
        P_Lattice = [1:1:size(p_Lattice,2);p_Lattice]';
        T_Lattice = [1:1:size(t_Lattice,2);t_Lattice]';
        File_ID = fopen('P_Lattice.txt','w+');
        for ii = 1:1:size(P_Lattice,1)
            fprintf(File_ID,'%d, %e, %e\n',P_Lattice(ii,1), ...
                                           P_Lattice(ii,2), ...
                                           P_Lattice(ii,3));
        end
        fclose(File_ID);
        File_ID = fopen('T_Lattice.txt','w+');
        for ii = 1:1:size(T_Lattice,1)
            fprintf(File_ID,'%d, %d, %d, %d\n',T_Lattice(ii,1), ...
                                               T_Lattice(ii,2), ...
                                               T_Lattice(ii,3), ...
                                               T_Lattice(ii,4));
        end
        fclose(File_ID);
        File_ID = fopen('Bottom.txt','w+');
        fprintf(File_ID,'%d, ',Bottom);
        fclose(File_ID);
        File_ID = fopen('Top.txt','w+');
        fprintf(File_ID,'%d, ',Top);
        fclose(File_ID);
        Abaqus_Input ( Material_Property );
        dos('abaqus job=Lattice.inp');
        fprintf('  -- Waiting for simulation completion ...\n');
        while double(isfile('Lattice.log'))==0
            pause(Pause_Time);
            fprintf('    -- Simulation is still running ...\n');
        end
        Flag = 0;
        while Flag == 0
            File_ID = fopen('Lattice.log');
            C = textscan(File_ID,'%s');
            fclose(File_ID);
            if cellfun(@numel, C) ~= 0
                if char(C{1}{cellfun(@numel, C)}) == "errors"
                    Flag = 1;
                    fprintf ...
                        ('    -- Simulation cannot run with errors ...\n');
                elseif char(C{1}{cellfun(@numel, C)}) == "COMPLETED"
                    Flag = 1;
                    fprintf('    -- Simulation is done ...\n');
                else
                    fprintf('    -- Simulation is still running ...\n');
                    Flag = 0;
                    pause(Pause_Time);
                end
            end
        end
        Text = fileread('Lattice.dat');
        Position = strfind(Text,'N O D E   O U T P U T');
        Force = abs(str2num(Text(Position+200:Position+215)));
        Stress = Force/Material_Property(1)/(Box_Crop(2)-Box_Crop(1));
        Strain = 1/(0.75*(Box_Crop(4)-Box_Crop(3)));
        E_Macro = Stress/Strain;
        cd ..
        pause(Pause_Time);
        rmdir 'Simulation_Folder' s;
        Feature_Matrix(Data_Count,:) = Feature;
        Property_Vector(Data_Count,1) = E_Macro/10000E6;
    end
    if mod(i,100) == 0
        if i == 100
            mkdir('Feature_Property_Output');
        end
        cd('Feature_Property_Output');
        dlmwrite(strcat('Feature_Matrix_',num2str(i/100), ...
                        '.txt'),Feature_Matrix);
        dlmwrite(strcat('Property_Vector_',num2str(i/100), ...
                        '.txt'),Property_Vector);
        clear Feature_Matrix Property_Vector
        Data_Count = 0;
        cd ..
    end
end
% =========================================================================




%  ============================== Subroutine ==============================
function  [ Line ] = Polygon_Line ( Node_A , Node_B , Thickness )
%  ------------------------------------------------------------------------
%  Input:  Node_A -- Coordinates of node A.
%          Node_B -- Coordinates of node B.
%          Thickness -- Thickness of line.
%  Output: Line -- Polygon of line.
%  ------------------------------------------------------------------------
Vector_Original = Node_B - Node_A;
Vector_Perpendicular = [-Vector_Original(2);Vector_Original(1)];
Arrow_Perpendicular = Vector_Perpendicular/norm(Vector_Perpendicular);
Point_1 = Node_A + 1/2*Thickness*Arrow_Perpendicular;
Point_2 = Node_A - 1/2*Thickness*Arrow_Perpendicular;
Point_3 = Node_B - 1/2*Thickness*Arrow_Perpendicular;
Point_4 = Node_B + 1/2*Thickness*Arrow_Perpendicular;
Line_X = [Point_1(1),Point_2(1),Point_3(1),Point_4(1)];
Line_Y = [Point_1(2),Point_2(2),Point_3(2),Point_4(2)];
Line = polyshape(Line_X,Line_Y);
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Polygon_Vector_Macro ] = Polygon_Union ( Polygon_Vector )
%  ------------------------------------------------------------------------
%  Input:  Polygon_Vector -- Vector of polygons.
%  Output: Polygon_Vector_Macro -- Macro polygon after union.
%  ------------------------------------------------------------------------
Sqrt_Number = ceil(sqrt(length(Polygon_Vector)));
Polygon_Vector_Macro = [];
for i = 1:1:ceil(length(Polygon_Vector)/Sqrt_Number)
    if i ~= ceil(length(Polygon_Vector)/Sqrt_Number)
        Polygon_Vector_Bulk = union(Polygon_Vector([1:Sqrt_Number] + ...
                                                   (i-1)*Sqrt_Number));
        Polygon_Vector_Macro = [Polygon_Vector_Macro,Polygon_Vector_Bulk];
    else
        Polygon_Vector_Bulk = union(Polygon_Vector ...
                                    ([(i-1)*Sqrt_Number+1:end]));
        Polygon_Vector_Macro = [Polygon_Vector_Macro,Polygon_Vector_Bulk];
    end
end
Polygon_Vector_Macro = union(Polygon_Vector_Macro);
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