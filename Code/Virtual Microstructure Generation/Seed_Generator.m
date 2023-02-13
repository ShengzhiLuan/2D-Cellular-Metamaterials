%  Generation of Seed Pattern by Domain, Seed Number, Type and Variation
%  Johns Hopkins University
%  Shengzhi Luan
%  02.08.2023
% =========================================================================
function  [ Seed ] = Seed_Generator ...
                                   ( Box , Seed_Number , Type , Variation )
%  ------------------------------------------------------------------------
%  Input:  Box -- Dimension of box ( [ x_min , x_max , y_min , y_max ] ).
%          Seed_Number -- Number of original seed pattern.
%          Type -- Lattice type definition:
%                  Type(1) -- Parameter k of k-uniform tiling: 
%                             ( 0 for random pattern,
%                               1 or 2 for 1- or 2-uniform tiling. ).
%                  Type(2) -- Index of uniform tiling:
%                             ( 1 -- [ 3 3 3 3 3 3
%                                      4 4 4 4
%                                      6 6 6
%                                      3 3 3 4 4
%                                      3 3 4 3 4
%                                      3 3 3 3 6
%                                      3 6 3 6
%                                      3 4 6 4
%                                      3 12 12
%                                      4 6 12
%                                      4 8 8 ], 
%                               2 -- [ 3 3 3 3 3 3 ; 3 3 3 4 4 (1)
%                                      3 3 3 3 3 3 ; 3 3 3 4 4 (2)
%                                      3 3 3 3 3 3 ; 3 3 4 3 4
%                                      3 3 3 3 3 3 ; 3 3 3 3 6 (1)
%                                      3 3 3 3 3 3 ; 3 3 3 3 6 (2)
%                                      3 3 3 3 3 3 ; 3 3 6 6
%                                      3 3 3 3 3 3 ; 3 3 4 12
%                                      4 4 4 4 ; 3 3 3 4 4 (1)
%                                      4 4 4 4 ; 3 3 3 4 4 (2)
%                                      3 3 3 4 4 ; 3 4 6 4
%                                      3 3 3 4 4 ; 3 3 4 3 4 (1)
%                                      3 3 3 4 4 ; 3 3 4 3 4 (2)
%                                      3 3 4 3 4 ; 3 4 6 4
%                                      3 3 3 3 6 ; 3 3 6 6
%                                      3 3 6 6 ; 3 6 3 6
%                                      3 6 3 6 ; 3 4 4 6 (1)
%                                      3 6 3 6 ; 3 4 4 6 (2)
%                                      3 4 4 6 ; 3 4 6 4
%                                      3 4 6 4 ; 4 6 12
%                                      3 4 3 12 ; 3 12 12 ]. ).
%                  Type(3) -- Dual pattern flag ( 0 for original pattern,
%                                                 1 for dual pattern. ).
%          Variation -- Variation coefficients:
%                       ( 1st row is seed perturbation,
%                         2nd row is global weight distinctiveness,
%                         3rd row is individual weight variance. ).
%  Output: Seed -- Seed pattern ( 1st row is x-coordinate,
%                                 2nd row is y-coordinate,
%                                 3rd row is weight. ).
%  ------------------------------------------------------------------------
Lx = Box(2)-Box(1);
Ly = Box(4)-Box(3);
Overlap_Coefficient = 0.3;
%  ------------------------------------------------------------------------
L_Min = 0;
L_Max = sqrt(Lx^2+Ly^2)/2;
Seed_Number_Trial = 0;
%  ------------------------------------------------------------------------
if Type(1) ~= 0
    count = 0;
    while Seed_Number_Trial ~= Seed_Number
        count = count + 1;
        Length = 1/2*(L_Min+L_Max);
        [ Seed ] = Seed_Pattern ( Box , Length , Type );
        Seed_Number_Trial = size(Seed,2);
        if Seed_Number_Trial > Seed_Number
            L_Min = Length;
        end
        if Seed_Number_Trial < Seed_Number
            L_Max = Length;
        end
        if count >= 100
            break
        end
    end
    if Type(3) == 1
        [ ~ , Seed , ~ , ~ , ~ ] = Seed_to_Tessellation ( Seed );
        Seed(3,:) = 1/2*Length;
    end
end
%  ------------------------------------------------------------------------
if Type(1) == 0
    L_Mean = sqrt(Lx*Ly/Seed_Number/pi);
    L_Var = Variation(3);
    Weight_Vector = normrnd(1,L_Var,[1,Seed_Number]);
    Weight_Vector = abs(sort(Weight_Vector,'Descend')*L_Mean);
    Overlap_Flag = zeros(1,length(Weight_Vector));
    Delete = [];
    for i = 1:1:length(Weight_Vector)
        count = 0;
        while Overlap_Flag(i) == 0
            count = count + 1;
            Overlap_Flag(i) = 1;
            Center(:,i) = [(Box(2)-Box(1))*rand+Box(1)
                           (Box(4)-Box(3))*rand+Box(3)];
            for j = 1:1:i-1
                if norm(Center(:,i)-Center(:,j)) + ...
                   Overlap_Coefficient*Weight_Vector(j) < ...
                   sum(Weight_Vector([i,j]))
                    Overlap_Flag(i) = 0;
                    break
                end
            end
            if count >= 100
                Delete(i) = 1;
                break
            end
        end
    end
    Center(:,find(Delete==1)) = [];
    Weight_Vector(:,find(Delete==1)) = [];
    Seed = [Center;Weight_Vector];
end
%  ------------------------------------------------------------------------
if Type(1) ~= 0
    [ Seed ] = Seed_Variation ( Seed , Variation );
end
%  ------------------------------------------------------------------------
end
% =========================================================================




%  ============================== Subroutine ==============================
function  [ Seed ] = Seed_Pattern ( Box , Length , Type )
%  ------------------------------------------------------------------------
%  Input:  Box -- Dimension of box ( [ x_min , x_max , y_min , y_max ] ).
%          Length -- Length of basic seed radius.
%          Type -- Lattice type definition:
%                  Type(1) -- Parameter k of k-uniform tiling: 
%                             ( 0 for random pattern,
%                               1 or 2 for 1- or 2-uniform tiling. ).
%                  Type(2) -- Index of uniform tiling:
%                             ( 1 -- [ 3 3 3 3 3 3
%                                      4 4 4 4
%                                      6 6 6
%                                      3 3 3 4 4
%                                      3 3 4 3 4
%                                      3 3 3 3 6
%                                      3 6 3 6
%                                      3 4 6 4
%                                      3 12 12
%                                      4 6 12
%                                      4 8 8 ], 
%                               2 -- [ 3 3 3 3 3 3 ; 3 3 3 4 4 (1)
%                                      3 3 3 3 3 3 ; 3 3 3 4 4 (2)
%                                      3 3 3 3 3 3 ; 3 3 4 3 4
%                                      3 3 3 3 3 3 ; 3 3 3 3 6 (1)
%                                      3 3 3 3 3 3 ; 3 3 3 3 6 (2)
%                                      3 3 3 3 3 3 ; 3 3 6 6
%                                      3 3 3 3 3 3 ; 3 3 4 12
%                                      4 4 4 4 ; 3 3 3 4 4 (1)
%                                      4 4 4 4 ; 3 3 3 4 4 (2)
%                                      3 3 3 4 4 ; 3 4 6 4
%                                      3 3 3 4 4 ; 3 3 4 3 4 (1)
%                                      3 3 3 4 4 ; 3 3 4 3 4 (2)
%                                      3 3 4 3 4 ; 3 4 6 4
%                                      3 3 3 3 6 ; 3 3 6 6
%                                      3 3 6 6 ; 3 6 3 6
%                                      3 6 3 6 ; 3 4 4 6 (1)
%                                      3 6 3 6 ; 3 4 4 6 (2)
%                                      3 4 4 6 ; 3 4 6 4
%                                      3 4 6 4 ; 4 6 12
%                                      3 4 3 12 ; 3 12 12 ]. ).
%                  Type(3) -- Dual pattern flag ( 0 for original pattern,
%                                                 1 for dual pattern. ).
%  Output: Seed -- Seed original pattern ( 1st row is x-coordinate,
%                                          2nd row is y-coordinate,
%                                          3rd row is weight. ).
%  ------------------------------------------------------------------------
x_min = Box(1);
x_max = Box(2);
y_min = Box(3);
y_max = Box(4);
x_range = x_max - x_min;
y_range = y_max - y_min;
%  ------------------------------------------------------------------------
if Type(1) == 1
    %  --------------------------------------------------------------------
    if Type(2) == 1
        count = 0;
        for j = 1:1:1E10
            y = (floor((j-1)/2)*sqrt(3)/2 + mod(j-1,2)*sqrt(3)/6)*Length;
            w = sqrt(3)/6*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,4) == 1 || mod(j,4) == 0
                    x = (i-1)*Length;
                end
                if mod(j,4) == 2 || mod(j,4) == 3
                    x = (i-1/2)*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 2
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*Length;
            w = 1/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                x = (i-1)*Length;
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 3
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*sqrt(3)/2*Length;
            w = sqrt(3)/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*3*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*3*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 4
        count = 0;
        for j = 1:1:1E10
            if mod(j,3) == 1
                y = floor((j-1)/3)*(1+sqrt(3)/2)*Length;
                w = 1/2*Length;
            end
            if mod(j,3) == 2
                y = (floor((j-1)/3)*(1+sqrt(3)/2) + 1/2+sqrt(3)/6)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,3) == 0
                y = (floor((j-1)/3)*(1+sqrt(3)/2) + 1/2+sqrt(3)/3)*Length;
                w = sqrt(3)/6*Length;
            end
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,6) == 1 || mod(j,6) == 2 || mod(j,6) == 0
                    x = (i-1)*Length;
                end
                if mod(j,6) == 3 || mod(j,6) == 4 || mod(j,6) == 5
                    x = (i-1/2)*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 5
        count = 0;
        for j = 1:1:1E10
            if mod(j,4) == 1
                y = floor((j-1)/4)*(1+sqrt(3))/2*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,4) == 2
                y = (floor((j-1)/4)*(1+sqrt(3))/2 + sqrt(3)/6)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,4) == 3
                y = (floor((j-1)/4)*(1+sqrt(3))/2 + (1+sqrt(3))/4)*Length;
                w = 1/2*Length;
            end
            if mod(j,4) == 0
                y = (floor((j-1)/4)*(1+sqrt(3))/2 + 1/2+sqrt(3)/3)*Length;
                w = sqrt(3)/6*Length;
            end
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,8) == 1
                    x = (floor((i-1)/2)*(1+sqrt(3)) + ...
                         mod(i-1,2)*sqrt(3)/3)*Length;
                end
                if mod(j,8) == 2 || mod(j,8) == 0
                    x = ((i-1)*(1+sqrt(3)) + 1/2+sqrt(3)*2/3)*Length;
                end
                if mod(j,8) == 3 || mod(j,8) == 7
                    x = ((i-1)*(1+sqrt(3))/2 + 1/4+sqrt(3)*5/12)*Length;
                end
                if mod(j,8) == 4 || mod(j,8) == 6
                    x = ((i-1)*(1+sqrt(3)) + sqrt(3)/6)*Length;
                end
                if mod(j,8) == 5
                    x = (floor((i-1)/2)*(1+sqrt(3)) + ...
                         mod(i-1,2)*sqrt(3)/3 + ...
                         1/2*(1+sqrt(3)))*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 6
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*sqrt(3)/2*Length;
            w = sqrt(3)/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                x = ((i-1)*7 + (j-1)*5/2 - floor((j-1)*5/2/7)*7)*Length;
                if x > x_range
                    break
                end
                count = count + 1;
                S1(:,count) = [x;y;w];
            end
        end
        CS = [];
        for i = 1:1:size(S1,2)
            C = S1(1:2,i);
            CX = [ 3 , 2 , 2 , 1 , 0 , -1 , -2 , -2 , -3 , ...
                  -3 , -2 , -2 , -1 , 0 , 1 , 2 , 2 , 3 ] * 1/2*Length;
            CY = [ 1 , 2 , 4 , 5 , 4 , 5 , 4 , 2 , 1 , ...
                  -1 , -2 , -4 , -5 , -4 , -5 , -4 , -2 , -1 ] * ...
                   sqrt(3)/6*Length;
            CW = ones(1,18)*sqrt(3)/6*Length;
            CS(1:2,(i-1)*18+1:(i-1)*18+18) = C + [CX;CY];
            CS(3,(i-1)*18+1:(i-1)*18+18) = CW;
        end
        S = [S1,CS];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 7
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*sqrt(3)/3*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*2*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*2*Length;
                end
                if mod(j,3) == 1 || mod(j,3) == 0
                    w = sqrt(3)/6*Length;
                end
                if mod(j,3) == 2
                    w = sqrt(3)/2*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 8
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*(1+sqrt(3))/2*Length;
            w = sqrt(3)/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*(3+sqrt(3))*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*(3+sqrt(3))*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S1(:,count) = [x;y;w];
            end
        end
        CS = [];
        for i = 1:1:size(S1,2)
            C = S1(1:2,i);
            CX = [ 1/3 , 1/4 , 1/6 , 0 , -1/6 , -1/4 , ...
                  -1/3 , -1/4 , -1/6 , 0 , 1/6 , 1/4 ] * ...
                   (3+sqrt(3))*Length;
            CY = [ 0 , 1/4 , 1/2 , 1/2 , 1/2 , 1/4 , ...
                   0 , -1/4 , -1/2 , -1/2 , -1/2 , -1/4 ] * ...
                   (1+sqrt(3))*Length;
            CW = [ sqrt(3)/6 , 1/2 , sqrt(3)/6 , 1/2 , ...
                   sqrt(3)/6 , 1/2 , sqrt(3)/6 , 1/2 , ...
                   sqrt(3)/6 , 1/2 , sqrt(3)/6 , 1/2 ] * Length;
            CS(1:2,(i-1)*12+1:(i-1)*12+12) = C + [CX;CY];
            CS(3,(i-1)*12+1:(i-1)*12+12) = CW;
        end
        S = [S1,CS];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 9
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*(1/2+sqrt(3)/3)*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*(2+sqrt(3))*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*(2+sqrt(3))*Length;
                end
                if mod(j,3) == 1 || mod(j,3) == 0
                    w = sqrt(3)/6*Length;
                end
                if mod(j,3) == 2
                    w = (2+sqrt(3))/2*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 10
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*(3+sqrt(3))/2*Length;
            w = (2+sqrt(3))/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*3*(1+sqrt(3))*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*3*(1+sqrt(3))*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S1(:,count) = [x;y;w];
            end
        end
        CS1 = [];
        CS2 = [];
        for i = 1:1:size(S1,2)
            C = S1(1:2,i);
            R1 = (1+sqrt(3))*Length;
            Theta1 = [0:60:300]/180*pi;
            CS1(1:2,(i-1)*6+1:(i-1)*6+6) = C+R1*[cos(Theta1);sin(Theta1)];
            CS1(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*sqrt(3)/2*Length;
            R2 = (3+sqrt(3))/2*Length;
            Theta2 = [30:60:330]/180*pi;
            CS2(1:2,(i-1)*6+1:(i-1)*6+6) = C+R2*[cos(Theta2);sin(Theta2)];
            CS2(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*1/2*Length;
        end
        S = [S1,CS1,CS2];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 11
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*(1+sqrt(2)/2)*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                x = (i-1)*(1+sqrt(2)/2)*Length;
                if mod(j,2) + mod(i,2) == 1
                    w = 1/2*Length;
                else
                    w = (1+sqrt(2))/2*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
end
%  ------------------------------------------------------------------------
if Type(1) == 2
    %  --------------------------------------------------------------------
    if Type(2) == 1
        count = 0;
        for j = 1:1:1E10
            if mod(j,5) == 1
                y = floor((j-1)/5)*(1+sqrt(3))*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,5) == 2
                y = (floor((j-1)/5)*(1+sqrt(3)) + sqrt(3)/6)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,5) == 3
                y = (floor((j-1)/5)*(1+sqrt(3)) + sqrt(3)/2)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,5) == 4
                y = (floor((j-1)/5)*(1+sqrt(3)) + sqrt(3)*2/3)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,5) == 0
                y = (floor((j-1)/5)*(1+sqrt(3)) + 1/2+sqrt(3)*5/6)*Length;
                w = 1/2*Length;
            end
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,5) == 1 || mod(j,5) == 4 || mod(j,5) == 0
                    x = (i-1)*Length;
                end
                if mod(j,5) == 2 || mod(j,5) == 3
                    x = (i-1/2)*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 2
        count = 0;
        for j = 1:1:1E10
            if mod(j,7) == 1
                y = floor((j-1)/7)*(1+sqrt(3)*3/2)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,7) == 2
                y = (floor((j-1)/7)*(1+sqrt(3)*3/2) + sqrt(3)/6)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,7) == 3
                y = (floor((j-1)/7)*(1+sqrt(3)*3/2) + sqrt(3)/2)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,7) == 4
                y = (floor((j-1)/7)*(1+sqrt(3)*3/2) + sqrt(3)*2/3)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,7) == 5
                y = (floor((j-1)/7)*(1+sqrt(3)*3/2) + sqrt(3))*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,7) == 6
                y = (floor((j-1)/7)*(1+sqrt(3)*3/2) + sqrt(3)*7/6)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,7) == 0
                y = (floor((j-1)/7)*(1+sqrt(3)*3/2) + ...
                     1/2+sqrt(3)*4/3)*Length;
                w = 1/2*Length;
            end
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,14) == 1 || mod(j,14) == 4 || ...
                   mod(j,14) == 5 || mod(j,14) == 9 || ...
                   mod(j,14) == 10 || mod(j,14) == 13 || mod(j,14) == 0
                    x = (i-1)*Length;
                end
                if mod(j,14) == 2 || mod(j,14) == 3 || ...
                   mod(j,14) == 6 || mod(j,14) == 7 || ...
                   mod(j,14) == 8 || mod(j,14) == 11 || mod(j,14) == 12
                    x = (i-1/2)*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 3
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*(1+sqrt(3))/2*Length;
            w = 0;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*(3+sqrt(3))*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*(3+sqrt(3))*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S0(:,count) = [x;y;w];
            end
        end
        CS1 = [];
        CS2 = [];
        for i = 1:1:size(S0,2)
            C = S0(1:2,i);
            CX = [ 1/3 , 1/4 , 1/6 , 0 , -1/6 , -1/4 , ...
                  -1/3 , -1/4 , -1/6 , 0 , 1/6 , 1/4 ] * ...
                   (3+sqrt(3))*Length;
            CY = [ 0 , 1/4 , 1/2 , 1/2 , 1/2 , 1/4 , ...
                   0 , -1/4 , -1/2 , -1/2 , -1/2 , -1/4 ] * ...
                   (1+sqrt(3))*Length;
            CW = [ sqrt(3)/6 , 1/2 , sqrt(3)/6 , 1/2 , ...
                   sqrt(3)/6 , 1/2 , sqrt(3)/6 , 1/2 , ...
                   sqrt(3)/6 , 1/2 , sqrt(3)/6 , 1/2 ] * Length;
            CS1(1:2,(i-1)*12+1:(i-1)*12+12) = C + [CX;CY];
            CS1(3,(i-1)*12+1:(i-1)*12+12) = CW;
            R2 = sqrt(3)/3*Length;
            Theta2 = [30:60:330]/180*pi;
            CS2(1:2,(i-1)*6+1:(i-1)*6+6) = C+R2*[cos(Theta2);sin(Theta2)];
            CS2(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*sqrt(3)/6*Length;
        end
        S = [CS1,CS2];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 4
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*sqrt(3)*3/2*Length;
            w = sqrt(3)/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*3*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*3*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S1(:,count) = [x;y;w];
            end
        end
        CS = [];
        for i = 1:1:size(S1,2)
            C = S1(1:2,i);
            CX = [ 3 , 2 , 2 , 1 , 0 , -1 , -2 , -2 , -3 , ...
                  -3 , -2 , -2 , -1 , 0 , 1 , 2 , 2 , 3 ] * 1/2*Length;
            CY = [ 1 , 2 , 4 , 5 , 4 , 5 , 4 , 2 , 1 , ...
                  -1 , -2 , -4 , -5 , -4 , -5 , -4 , -2 , -1 ] * ...
                   sqrt(3)/6*Length;
            CW = ones(1,18)*sqrt(3)/6*Length;
            CS(1:2,(i-1)*18+1:(i-1)*18+18) = C + [CX;CY];
            CS(3,(i-1)*18+1:(i-1)*18+18) = CW;
        end
        S = [S1,CS];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 5
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*sqrt(3)/2*Length;
            w = sqrt(3)/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                x = ((i-1)*13 + (j-1)*7/2 - floor((j-1)*7/2/13)*13)*Length;
                if x > x_range
                    break
                end
                count = count + 1;
                S1(:,count) = [x;y;w];
            end
        end
        CS = [];
        for i = 1:1:size(S1,2)
            C = S1(1:2,i);
            CX = [ 5 , 4 , 3 , 2 , 4 , 3 , 2 , 1 , ...
                   0 , 3 , 2 , 1 , 0 , -1, -2 , -3 , ...
                  -1 , -2 , -3 , -4 , -2 , -3 , -4 , -5 , ...
                  -5 , -4 , -3 , -2 , -4 , -3 , -2 , -1 , ...
                   0 , -3, -2 , -1 , 0 , 1 , 2 , 3 , ...
                   1 , 2 , 3 , 4 , 2 , 3 , 4 , 5 ] * 1/2*Length;
            CY = [ 1 , 2 , 1 , 2 , 4 , 5 , 4 , 5 , ...
                   4 , 7 , 8 , 7 , 8 , 7 , 8 , 7 , ...
                   5 , 4 , 5 , 4 , 2 , 1 , 2 , 1 , ...
                  -1 , -2 , -1 , -2 , -4 , -5 , -4 , -5 , ...
                  -4 , -7 , -8 , -7 , -8 , -7 , -8 , -7 , ...
                  -5 , -4 , -5 , -4 , -2 , -1 , -2 , -1 ] * ...
                  sqrt(3)/6*Length;
            CW = ones(1,48)*sqrt(3)/6*Length;
            CS(1:2,(i-1)*48+1:(i-1)*48+48) = C + [CX;CY];
            CS(3,(i-1)*48+1:(i-1)*48+48) = CW;
        end
        S = [S1,CS];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 6
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*sqrt(3)*3/2*Length;
            w = 0;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*3*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*3*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S0(:,count) = [x;y;w];
            end
        end
        CS1 = [];
        CS2 = [];
        for i = 1:1:size(S0,2)
            C = S0(1:2,i);
            R1 = sqrt(3)/3*Length;
            Theta1 = [30:60:330]/180*pi;
            CS1(1:2,(i-1)*6+1:(i-1)*6+6) = C+R1*[cos(Theta1);sin(Theta1)];
            CS1(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*sqrt(3)/6*Length;
            R2 = sqrt(3)*Length;
            Theta2 = [30:60:330]/180*pi;
            CS2(1:2,(i-1)*6+1:(i-1)*6+6) = C+R2*[cos(Theta2);sin(Theta2)];
            CS2(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*sqrt(3)/2*Length;
        end
        S = [CS1,CS2];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 7
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*(3+sqrt(3))/2*Length;
            w = (2+sqrt(3))/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*3*(1+sqrt(3))*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*3*(1+sqrt(3))*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S1(:,count) = [x;y;w];
            end
        end
        CS11 = [];
        CS2 = [];
        for i = 1:1:size(S1,2)
            C = S1(1:2,i);
            R1 = (1+sqrt(3))*Length;
            Theta1 = [0:60:300]/180*pi;
            CS1(1:2,(i-1)*6+1:(i-1)*6+6) = C+R1*[cos(Theta1);sin(Theta1)];
            CS1(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*0;
            R2 = (3+sqrt(3))/2*Length;
            Theta2 = [30:60:330]/180*pi;
            CS2(1:2,(i-1)*6+1:(i-1)*6+6) = C+R2*[cos(Theta2);sin(Theta2)];
            CS2(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*1/2*Length;
        end
        for i = 1:1:size(CS1,2)
            C1 = CS1(1:2,i);
            R11 = sqrt(3)/3*Length;
            Theta11 = [0:60:300]/180*pi;
            CS11(1:2,(i-1)*6+1:(i-1)*6+6) = C1 + R11*[cos(Theta11)
                                                      sin(Theta11)];
            CS11(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*sqrt(3)/6*Length;
        end
        S = [S1,CS11,CS2];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 8
        count = 0;
        for j = 1:1:1E10
            if mod(j,4) == 1
                y = floor((j-1)/4)*(2+sqrt(3)/2)*Length;
                w = 1/2*Length;
            end
            if mod(j,4) == 2
                y = (floor((j-1)/4)*(2+sqrt(3)/2) + 1)*Length;
                w = 1/2*Length;
            end
            if mod(j,4) == 3
                y = (floor((j-1)/4)*(2+sqrt(3)/2) + 3/2+sqrt(3)/6)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,4) == 0
                y = (floor((j-1)/4)*(2+sqrt(3)/2) + 3/2+sqrt(3)/3)*Length;
                w = sqrt(3)/6*Length;
            end
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,8) <= 3.5
                    x = (i-1)*Length;
                end
                if mod(j,8) >= 3.5
                    x = (i-1/2)*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 9
        count = 0;
        for j = 1:1:1E10
            if mod(j,5) == 1
                y = floor((j-1)/5)*(3+sqrt(3)/2)*Length;
                w = 1/2*Length;
            end
            if mod(j,5) == 2
                y = (floor((j-1)/5)*(3+sqrt(3)/2) + 1)*Length;
                w = 1/2*Length;
            end
            if mod(j,5) == 3
                y = (floor((j-1)/5)*(3+sqrt(3)/2) + 2)*Length;
                w = 1/2*Length;
            end
            if mod(j,5) == 4
                y = (floor((j-1)/5)*(3+sqrt(3)/2) + 5/2+sqrt(3)/6)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,5) == 0
                y = (floor((j-1)/5)*(3+sqrt(3)/2) + 5/2+sqrt(3)/3)*Length;
                w = sqrt(3)/6*Length;
            end
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,10) <= 4.5
                    x = (i-1)*Length;
                end
                if mod(j,10) >= 4.5
                    x = (i-1/2)*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 10
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*(1+sqrt(3)/2)*Length;
            w = sqrt(3)/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*(3+2*sqrt(3))*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*(3+2*sqrt(3))*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S1(:,count) = [x;y;w];
            end
        end
        CS1 = [];
        CS2 = [];
        CS3 = [];
        CS31 = [];
        for i = 1:1:size(S1,2)
            C = S1(1:2,i);
            R1 = (1+sqrt(3))/2*Length;
            Theta1 = [30:60:330]/180*pi;
            CS1(1:2,(i-1)*6+1:(i-1)*6+6) = C+R1*[cos(Theta1);sin(Theta1)];
            CS1(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*1/2*Length;
            R2 = (3+sqrt(3))/2*Length;
            Theta2 = [30:60:330]/180*pi;
            CS2(1:2,(i-1)*6+1:(i-1)*6+6) = C+R2*[cos(Theta2);sin(Theta2)];
            CS2(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*1/2*Length;
            R3 = (1+sqrt(3)*2/3)*Length;
            Theta3 = [0:60:300]/180*pi;
            CS3(1:2,(i-1)*6+1:(i-1)*6+6) = C+R3*[cos(Theta3);sin(Theta3)];
            CS3(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*sqrt(3)/6*Length;
        end
        for i = 1:1:size(CS3,2)
            C3 = CS3(1:2,i);
            R31 = sqrt(3)/3*Length;
            Theta31 = ([60:120:300]+mod(i-1,6)*60)/180*pi;
            CS31(1:2,(i-1)*3+1:(i-1)*3+3) = C3+R31*[cos(Theta31)
                                                   sin(Theta31)];
            CS31(3,(i-1)*3+1:(i-1)*3+3) = ones(1,3)*sqrt(3)/6*Length;
        end
        S = [S1,CS1,CS2,CS3,CS31];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 11
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*(3+sqrt(3))/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                x = (i-1)*(3+sqrt(3))/2*Length;
                if x > x_range
                    break
                end
                count = count + 1;
                if mod(j,2) == 1
                    if mod(i,2) == 1
                        w = 1;
                    else
                        w = 2;
                    end
                end
                if mod(j,2) == 0
                    if mod(i,2) == 1
                        w = 2;
                    else
                        w = 1;
                    end
                end
                S0(:,count) = [x;y;w];
            end
        end
        CS = [];
        for i = 1:1:size(S0,2)
            C = S0(1:2,i);
            if S0(3,i) == 1
                CX = [ 1/2 , 1/2 , 0 , -1/2 , -1/2 , -1/2 , 0 , 1/2 , ...
                       1+sqrt(3)/6 , 1+sqrt(3)/3 , (3+sqrt(3))/4 , ...
                       1/2+sqrt(3)/6 , 0 , -1/2-sqrt(3)/6 , ...
                      -(3+sqrt(3))/4 , -1-sqrt(3)/3 , -1-sqrt(3)/6 ...
                      -1-sqrt(3)/3 , -(3+sqrt(3))/4 , -1/2-sqrt(3)/6 , ...
                       0 , 1/2+sqrt(3)/6 , (3+sqrt(3))/4 , ...
                       1+sqrt(3)/3 ] * Length;
                CY = [ 0 , 1/2+sqrt(3)/6 , 1/2+sqrt(3)/3 , ...
                       1/2+sqrt(3)/6 , 0 , -1/2-sqrt(3)/6 , ...
                      -1/2-sqrt(3)/3 , -1/2-sqrt(3)/6 , 0 , 1/2 , ...
                       (3+sqrt(3))/4 , 1+sqrt(3)/2 , 1+sqrt(3)/2 , ...
                       1+sqrt(3)/2 , (3+sqrt(3))/4 , 1/2 , 0 , -1/2 , ...
                      -(3+sqrt(3))/4 , -1-sqrt(3)/2 , -1-sqrt(3)/2 , ...
                      -1-sqrt(3)/2 , -(3+sqrt(3))/4 , -1/2 ] * Length;
                CW = [ 1/2 , sqrt(3)/6 , sqrt(3)/6 , sqrt(3)/6 , ...
                       1/2 , sqrt(3)/6 , sqrt(3)/6 , sqrt(3)/6 , ...
                       sqrt(3)/6 , sqrt(3)/6 , 1/2 , sqrt(3)/6 , ...
                       1/2 , sqrt(3)/6 , 1/2 , sqrt(3)/6 , ...
                       sqrt(3)/6 , sqrt(3)/6 , 1/2 , sqrt(3)/6 , ...
                       1/2 , sqrt(3)/6 , 1/2 , sqrt(3)/6 ] * Length;
                CS(1:2,(i-1)*24+1:(i-1)*24+24) = C + [CX;CY];
                CS(3,(i-1)*24+1:(i-1)*24+24) = CW;
            end
            if S0(3,i) == 2
                CX = [ 0 , 1/2+sqrt(3)/6 , 1/2+sqrt(3)/3 , ...
                       1/2+sqrt(3)/6 , 0 , -1/2-sqrt(3)/6 , ...
                      -1/2-sqrt(3)/3 , -1/2-sqrt(3)/6 , 0 , 1/2 , ...
                       (3+sqrt(3))/4 , 1+sqrt(3)/2 , 1+sqrt(3)/2 , ...
                       1+sqrt(3)/2 , (3+sqrt(3))/4 , 1/2 , 0 , -1/2 , ...
                      -(3+sqrt(3))/4 , -1-sqrt(3)/2 , -1-sqrt(3)/2 , ...
                      -1-sqrt(3)/2 , -(3+sqrt(3))/4 , -1/2 ] * Length;
                CY = [ 1/2 , 1/2 , 0 , -1/2 , -1/2 , -1/2 , 0 , 1/2 , ...
                       1+sqrt(3)/6 , 1+sqrt(3)/3 , (3+sqrt(3))/4 , ...
                       1/2+sqrt(3)/6 , 0 , -1/2-sqrt(3)/6 , ...
                      -(3+sqrt(3))/4 , -1-sqrt(3)/3 , -1-sqrt(3)/6 ...
                      -1-sqrt(3)/3 , -(3+sqrt(3))/4 , -1/2-sqrt(3)/6 , ...
                       0 , 1/2+sqrt(3)/6 , (3+sqrt(3))/4 , ...
                       1+sqrt(3)/3 ] * Length;
                CW = [ 1/2 , sqrt(3)/6 , sqrt(3)/6 , sqrt(3)/6 , ...
                       1/2 , sqrt(3)/6 , sqrt(3)/6 , sqrt(3)/6 , ...
                       sqrt(3)/6 , sqrt(3)/6 , 1/2 , sqrt(3)/6 , ...
                       1/2 , sqrt(3)/6 , 1/2 , sqrt(3)/6 , ...
                       sqrt(3)/6 , sqrt(3)/6 , 1/2 , sqrt(3)/6 , ...
                       1/2 , sqrt(3)/6 , 1/2 , sqrt(3)/6 ] * Length;
                CS(1:2,(i-1)*24+1:(i-1)*24+24) = C + [CX;CY];
                CS(3,(i-1)*24+1:(i-1)*24+24) = CW;
            end
        end
        S = CS;
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 12
        S_Old = 1/2*[x_range;y_range;-2];
        Expand_ID = 1;
        Expand_Flag = 1;
        while Expand_Flag == 1
            count = 0;
            for i = 1:1:length(Expand_ID)
                if S_Old(3,Expand_ID(i)) == -1
                    count = count + 1;
                    X_Expand = [ 5/4+sqrt(3)/2 , -3/4-sqrt(3)/2 , ...
                                -5/4-sqrt(3)/2 , 3/4+sqrt(3)/2 ] * Length;
                    Y_Expand = [ sqrt(3)/4 , 1+sqrt(3)/4 , ...
                                -sqrt(3)/4 , -1-sqrt(3)/4 ] * Length;
                    W_Expand = [ -2 , -2 , -2 , -2 ];
                    S_Expand(1:2,(count-1)*4+1:(count-1)*4+4) = ...
                             S_Old(1:2,Expand_ID(i)) + [X_Expand;Y_Expand];
                    S_Expand(3,(count-1)*4+1:(count-1)*4+4) = W_Expand;
                end
                if S_Old(3,Expand_ID(i)) == -2
                    count = count + 1;
                    X_Expand = [ 5/4+sqrt(3)/2 , -3/4-sqrt(3)/2 , ...
                                -5/4-sqrt(3)/2 , 3/4+sqrt(3)/2 ] * Length;
                    Y_Expand = [ sqrt(3)/4 , 1+sqrt(3)/4 , ...
                                -sqrt(3)/4 , -1-sqrt(3)/4 ] * Length;
                    W_Expand = [ -1 , -1 , -1 , -1 ];
                    S_Expand(1:2,(count-1)*4+1:(count-1)*4+4) = ...
                             S_Old(1:2,Expand_ID(i)) + [X_Expand;Y_Expand];
                    S_Expand(3,(count-1)*4+1:(count-1)*4+4) = W_Expand;
                end
            end
            [~,Index,~] = uniquetol(S_Expand',1E-3,'Byrows',true);
            S_Expand = S_Expand(:,Index);
            count = 0;
            Delete_ID = [];
            for i = 1:1:size(S_Expand,2)
                if S_Expand(1,i) < 0 || S_Expand(1,i) > x_range || ...
                   S_Expand(2,i) < 0 || S_Expand(2,i) > y_range
                    count = count + 1;
                    Delete_ID(count) = i;
                end
                Distance_Vector = S_Expand(1:2,i) - S_Old(1:2,:);
                Distance = [];
                for j = 1:1:size(Distance_Vector,2)
                    Distance(j) = norm(Distance_Vector(:,j));
                end
                if min(Distance) <= 1E-5*Length
                    count = count + 1;
                    Delete_ID(count) = i;
                end
            end
            S_Expand(:,Delete_ID) = [];
            S_New = [S_Old,S_Expand];
            Expand_ID = [1:1:size(S_Expand,2)]+size(S_Old,2);
            S_Old = S_New;
            if length(Expand_ID) == 0
                Expand_Flag = 0;
            end
        end
        S0 = S_Old;
        CS = [];
        for i = 1:1:size(S0,2)
            C = S0(1:2,i);
            CX = [ 1/2 , 1/2 , 0 , -1/2 , -1/2 , -1/2 , 0 , 1/2 ]*Length;
            CY = [ 0 , 1/2+sqrt(3)/6 , 1/2+sqrt(3)/3 , 1/2+sqrt(3)/6 , ...
                   0 , -1/2-sqrt(3)/6 , -1/2-sqrt(3)/3 , -1/2-sqrt(3)/6 ...
                   ] * Length;
            CW = [ 1/2 , sqrt(3)/6 , sqrt(3)/6 , sqrt(3)/6 , ...
                   1/2 , sqrt(3)/6 , sqrt(3)/6 , sqrt(3)/6 ] * Length;
            Theta1 = 0/180*pi;
            Theta2 = -30/180*pi;
            if S0(3,i) == -1
                CS(1:2,(i-1)*8+1:(i-1)*8+8) = C + [CX;CY];
                CS(3,(i-1)*8+1:(i-1)*8+8) = CW;
            end
            if S0(3,i) == -2
                CS(1:2,(i-1)*8+1:(i-1)*8+8) = C + [CX*cos(Theta2) - ...
                                                   CY*sin(Theta2)
                                                   CX*sin(Theta2) + ...
                                                   CY*cos(Theta2)];
                CS(3,(i-1)*8+1:(i-1)*8+8) = CW;
            end
        end
        S = CS;
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 13
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*(1/2+sqrt(3)/3)*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*(2+sqrt(3))*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*(2+sqrt(3))*Length;
                end
                if mod(j,3) == 1 || mod(j,3) == 0
                    w = sqrt(3)/6*Length;
                end
                if mod(j,3) == 2
                    w = sqrt(3)/2*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S1(:,count) = [x;y;w];
            end
        end
        count = 0;
        CS1 = [];
        CS2 = [];
        for i = 1:1:size(S1,2)
            if S1(3,i) == sqrt(3)/2*Length
                count = count + 1;
                C = S1(1:2,i);
                R1 = (1+sqrt(3)/3)*Length;
                Theta1 = [0:60:300]/180*pi;
                CS1(1:2,(count-1)*6+1:(count-1)*6+6) = C+R1*[cos(Theta1)
                                                             sin(Theta1)];
                CS1(3,(count-1)*6+1:(count-1)*6+6) = ones(1,6) * ...
                                                     sqrt(3)/6*Length;
                R2 = (1+sqrt(3))/2*Length;
                Theta2 = [30:60:330]/180*pi;
                CS2(1:2,(count-1)*6+1:(count-1)*6+6) = C+R2*[cos(Theta2)
                                                             sin(Theta2)];
                CS2(3,(count-1)*6+1:(count-1)*6+6) = ones(1,6)*1/2*Length;
            end
        end
        S = [S1,CS1,CS2];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 14
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*5/2*Length;
            w = sqrt(3)/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*sqrt(3)*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*sqrt(3)*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S1(:,count) = [x;y;w];
            end
        end
        CS = [];
        for i = 1:1:size(S1,2)
            C = S1(1:2,i);
            CX = [ 2 , 1 , -1 , -2 , -2 , -1 , 1 , 2 ]*sqrt(3)/6*Length;
            CY = [ 2 , 3 , 3 , 2 , -2 , -3 , -3 , -2 ]*1/2*Length;
            CW = ones(1,8)*sqrt(3)/6*Length;
            CS(1:2,(i-1)*8+1:(i-1)*8+8) = C + [CX;CY];
            CS(3,(i-1)*8+1:(i-1)*8+8) = CW;
        end
        S = [S1,CS];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 15
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*sqrt(3)/3*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,3) == 1 || mod(j,3) == 0
                    x = (i-1)*2*Length;
                    w = sqrt(3)/6*Length;
                end
                if mod(j,3) == 2
                    x = (i-1/2)*2*Length;
                    w = sqrt(3)/2*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 16
        count = 0;
        for j = 1:1:1E10
            if mod(j,4) == 1
                y = floor((j-1)/4)*(1+sqrt(3))*Length;
            end
            if mod(j,4) == 2
                y = (floor((j-1)/4)*(1+sqrt(3)) + sqrt(3)/3)*Length;
            end
            if mod(j,4) == 3
                y = (floor((j-1)/4)*(1+sqrt(3)) + sqrt(3)*2/3)*Length;
            end
            if mod(j,4) == 0
                y = (floor((j-1)/4)*(1+sqrt(3)) + 1/2+sqrt(3)*5/6)*Length;
            end
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,4) == 1 || mod(j,4) == 3
                    x = (i-1)*2*Length;
                    w = sqrt(3)/6*Length;
                end
                if mod(j,4) == 2
                    x = (i-1/2)*2*Length;
                    w = sqrt(3)/2*Length;
                end
                if mod(j,4) == 0
                    x = (i-1)*Length;
                    w = 1/2*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 17
        count = 0;
        for j = 1:1:1E10
            if mod(j,4) == 1
                y = floor((j-1)/4)*(1+sqrt(3))*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,4) == 2
                y = (floor((j-1)/4)*(1+sqrt(3)) + sqrt(3)/3)*Length;
                w = sqrt(3)/2*Length;
            end
            if mod(j,4) == 3
                y = (floor((j-1)/4)*(1+sqrt(3)) + sqrt(3)*2/3)*Length;
                w = sqrt(3)/6*Length;
            end
            if mod(j,4) == 0
                y = (floor((j-1)/4)*(1+sqrt(3)) + 1/2+sqrt(3)*5/6)*Length;
                w = 1/2*Length;
            end
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,8) == 1 || mod(j,8) == 3 || mod(j,8) == 6
                    x = (i-1)*2*Length;
                end
                if mod(j,8) == 2 || mod(j,8) == 5 || mod(j,8) == 7
                    x = (i-1/2)*2*Length;
                end
                if mod(j,8) == 4 || mod(j,8) == 0
                    x = (i-1)*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
    if Type(2) == 18
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*(1+sqrt(3))*3/2*Length;
            w = sqrt(3)/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*(3+sqrt(3))*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*(3+sqrt(3))*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S1(:,count) = [x;y;w];
            end
        end
        CS1 = [];
        CS2 = [];
        CS3 = [];
        CS4 = [];
        for i = 1:1:size(S1,2)
            C = S1(1:2,i);
            R1 = (1+sqrt(3))/2*Length;
            Theta1 = [0:60:300]/180*pi;
            CS1(1:2,(i-1)*6+1:(i-1)*6+6) = C+R1*[cos(Theta1);sin(Theta1)];
            CS1(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*1/2*Length;
            R2 = (1+sqrt(3)/3)*Length;
            Theta2 = [30:60:330]/180*pi;
            CS2(1:2,(i-1)*6+1:(i-1)*6+6) = C+R2*[cos(Theta2);sin(Theta2)];
            CS2(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*sqrt(3)/6*Length;
            R3 = (3+sqrt(3))/2*Length;
            Theta3 = [0:60:300]/180*pi;
            CS3(1:2,(i-1)*6+1:(i-1)*6+6) = C+R3*[cos(Theta3);sin(Theta3)];
            CS3(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*1/2*Length;
            R4 = (1+sqrt(3))*Length;
            Theta4 = [30:60:330]/180*pi;
            CS4(1:2,(i-1)*6+1:(i-1)*6+6) = C+R4*[cos(Theta4);sin(Theta4)];
            CS4(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*sqrt(3)/2*Length;
        end
        S = [S1,CS1,CS2,CS3,CS4];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 19
        count = 0;
        for j = 1:1:1E10
            y = (j-1)*(3+sqrt(3))*Length;
            w = (2+sqrt(3))/2*Length;
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,2) == 1
                    x = (i-1)*2*(1+sqrt(3))*Length;
                end
                if mod(j,2) == 0
                    x = (i-1/2)*2*(1+sqrt(3))*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S1(:,count) = [x;y;w];
            end
        end
        CS1 = [];
        CS2 = [];
        CS3 = [];
        for i = 1:1:size(S1,2)
            C = S1(1:2,i);
            R1 = (1+sqrt(3))*Length;
            Theta1 = [0:60:300]/180*pi;
            CS1(1:2,(i-1)*6+1:(i-1)*6+6) = C+R1*[cos(Theta1);sin(Theta1)];
            CS1(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*sqrt(3)/2*Length;
            R2 = (3+sqrt(3))/2*Length;
            Theta2 = [30:60:330]/180*pi;
            CS2(1:2,(i-1)*6+1:(i-1)*6+6) = C+R2*[cos(Theta2);sin(Theta2)];
            CS2(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*1/2*Length;
            R3 = (2+sqrt(3)*2/3)*Length;
            Theta3 = [30:60:330]/180*pi;
            CS3(1:2,(i-1)*6+1:(i-1)*6+6) = C+R3*[cos(Theta3);sin(Theta3)];
            CS3(3,(i-1)*6+1:(i-1)*6+6) = ones(1,6)*sqrt(3)/6*Length;
        end
        S = [S1,CS1,CS2,CS3];
        [~,Index,~] = uniquetol(S(1:2,:)',1E-3,'Byrows',true);
        S = S(:,Index);
    end
    %  --------------------------------------------------------------------
    if Type(2) == 20
        count = 0;
        for j = 1:1:1E10
            if mod(j,4) == 1
                y = floor((j-1)/4)*(2+sqrt(3))*Length;
            end
            if mod(j,4) == 2
                y = (floor((j-1)/4)*(2+sqrt(3)) + 1/2+sqrt(3)/6)*Length;
            end
            if mod(j,4) == 3
                y = (floor((j-1)/4)*(2+sqrt(3)) + 1+sqrt(3)/2)*Length;
            end
            if mod(j,4) == 0
                y = (floor((j-1)/4)*(2+sqrt(3)) + 3/2+sqrt(3)*5/6)*Length;
            end
            if y > y_range
                break
            end
            for i = 1:1:1E10
                if mod(j,4) == 1
                    if mod(i,3) == 1
                        x = floor((i-1)/3)*(2+sqrt(3))*Length;
                        w = 1/2*Length;
                    end
                    if mod(i,3) == 2
                        x = (floor((i-1)/3)*(2+sqrt(3)) + ...
                             1/2+sqrt(3)/6)*Length;
                        w = sqrt(3)/6*Length;
                    end
                    if mod(i,3) == 0
                        x = (floor((i-1)/3)*(2+sqrt(3)) + ...
                             3/2+sqrt(3)*5/6)*Length;
                        w = sqrt(3)/6*Length;
                    end
                end
                if mod(j,4) == 2 || mod(j,4) == 0
                    x = (i-1)*(2+sqrt(3))*Length;
                    w = sqrt(3)/6*Length;
                end
                if mod(j,4) == 3
                    x = (i-1/2)*(2+sqrt(3))*Length;
                    w = (2+sqrt(3))/2*Length;
                end
                if x > x_range
                    break
                end
                count = count + 1;
                S(:,count) = [x;y;w];
            end
        end
    end
    %  --------------------------------------------------------------------
end
%  ------------------------------------------------------------------------
Seed = S + [x_min;y_min;0];
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Seed_Output ] = Seed_Variation ( Seed_Input , Variation )
%  ------------------------------------------------------------------------
%  Input:  Seed_Input -- Seed original pattern ( 1st row is x-coordinate,
%                                                2nd row is y-coordinate,
%                                                3rd row is weight. ).
%          Variation -- Variation coefficients:
%                       ( 1st row is seed perturbation,
%                         2nd row is global weight distinctiveness,
%                         3rd row is individual weight variance. ).
%  Output: Seed_Output -- Seed pattern ( 1st row is x-coordinate,
%                                        2nd row is y-coordinate,
%                                        3rd row is weight. ).
%  ------------------------------------------------------------------------
Seed = Seed_Input;
Seed_Coefficient = Variation(1);
GWD_Coefficient = Variation(2);
IWV_Coefficient = Variation(3);
%  ------------------------------------------------------------------------
count = 0;
for i = 1:1:9
    for j = i+1:1:10
        count = count + 1;
        Seed_Distance(count) = norm(Seed(1:2,i)-Seed(1:2,j));
    end
end
Length = min(Seed_Distance);
%  ------------------------------------------------------------------------
R_Coefficient = rand(1,size(Seed,2))*Seed_Coefficient*Length;
T_Coefficient = rand(1,size(Seed,2))*2*pi;
x_Coefficient = R_Coefficient.*cos(T_Coefficient);
y_Coefficient = R_Coefficient.*sin(T_Coefficient);
Seed(1:2,:) = Seed(1:2,:) + [x_Coefficient;y_Coefficient];
%  ------------------------------------------------------------------------
Seed(3,:) = Seed(3,:) + GWD_Coefficient*(mean(Seed(3,:))-Seed(3,:));
Seed(3,:) = Seed(3,:).*normrnd(1,IWV_Coefficient,[1,size(Seed,2)]);
%  ------------------------------------------------------------------------
Seed_Output = Seed;
%  ------------------------------------------------------------------------
end
%  ========================================================================