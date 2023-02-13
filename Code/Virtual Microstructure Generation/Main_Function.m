%  Generating Dataset of 2D Cellular Metamaterials
%  Johns Hopkins University
%  Shengzhi Luan
%  02.08.2023
% =========================================================================
Box = [0,1,0,1];
Min_Seed_Number = 1000;
Max_Seed_Number = 2500;
Lattice_Type_Matrix = [ 1 , 1 , 0
                        1 , 1 , 1 
                        1 , 2 , 0
                        1 , 2 , 1
                        1 , 3 , 0
                        1 , 3 , 1
                        1 , 4 , 0
                        1 , 4 , 1
                        1 , 5 , 0
                        1 , 5 , 1
                        1 , 6 , 0
                        1 , 6 , 1 
                        1 , 7 , 0
                        1 , 7 , 1
                        1 , 8 , 0
                        1 , 8 , 1
                        1 , 9 , 0
                        1 , 9 , 1
                        1 , 10 , 0
                        1 , 10 , 1
                        1 , 11 , 0
                        1 , 11 , 1
                        2 , 1 , 0
                        2 , 1 , 1 
                        2 , 2 , 0
                        2 , 2 , 1
                        2 , 3 , 0
                        2 , 3 , 1
                        2 , 4 , 0
                        2 , 4 , 1
                        2 , 5 , 0
                        2 , 5 , 1
                        2 , 6 , 0
                        2 , 6 , 1 
                        2 , 7 , 0
                        2 , 7 , 1
                        2 , 8 , 0
                        2 , 8 , 1
                        2 , 9 , 0
                        2 , 9 , 1
                        2 , 10 , 0
                        2 , 10 , 1
                        2 , 11 , 0
                        2 , 11 , 1
                        2 , 12 , 0
                        2 , 12 , 1
                        2 , 13 , 0
                        2 , 13 , 1 
                        2 , 14 , 0
                        2 , 14 , 1
                        2 , 15 , 0
                        2 , 15 , 1
                        2 , 16 , 0
                        2 , 16 , 1
                        2 , 17 , 0
                        2 , 17 , 1
                        2 , 18 , 0
                        2 , 18 , 1 
                        2 , 19 , 0
                        2 , 19 , 1
                        2 , 20 , 0
                        2 , 20 , 1
                        0 , 0 , 0 ];
Variation_Matrix = [ 0.20 , 0.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.40 , 0.00 , 0.40
                        0.40 , 0.00 , 0.40
                        0.30 , 0.00 , 0.30
                        0.20 , 0.00 , 0.20
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.30 , 1.00 , 0.30
                        0.30 , 0.00 , 0.30
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.30 , 1.00 , 0.30
                        0.30 , 0.00 , 0.30
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.40 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.40 , 1.00 , 0.40
                        0.30 , 0.00 , 0.30
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.30 , 1.00 , 0.30
                        0.30 , 0.00 , 0.30
                        0.30 , 1.00 , 0.30
                        0.20 , 0.00 , 0.30
                        0.30 , 1.00 , 0.30
                        0.20 , 0.00 , 0.30
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.30 , 1.00 , 0.30
                        0.30 , 0.00 , 0.30
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.40 , 1.00 , 0.30
                        0.30 , 0.00 , 0.30
                        0.20 , 1.00 , 0.30
                        0.30 , 0.00 , 0.30
                        0.30 , 1.00 , 0.30
                        0.30 , 0.00 , 0.30
                        0.30 , 1.00 , 0.30
                        0.30 , 0.00 , 0.30
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.20 , 1.00 , 0.20
                        0.30 , 0.00 , 0.30
                        0.30 , 1.00 , 0.30
                        0.30 , 0.00 , 0.30
                        0.40 , 1.00 , 0.20
                        0.20 , 0.00 , 0.30
                        0.40 , 1.00 , 0.15
                        0.30 , 0.00 , 0.30
                        0.00 , 0.00 , 0.30];
Sample_Number = [ 25 , 25 , 25 , 25 , 25 , 25 , 27 , 25 , 27 , 25 , ...
                  27 , 25 , 27 , 25 , 27 , 25 , 27 , 25 , 27 , 25 , ...
                  27 , 25 , 27 , 25 , 27 , 25 , 27 , 25 , 27 , 25 , ...
                  27 , 25 , 27 , 25 , 27 , 25 , 27 , 25 , 27 , 25 , ...
                  27 , 25 , 27 , 25 , 27 , 25 , 27 , 25 , 27 , 25 , ...
                  27 , 25 , 27 , 25 , 27 , 25 , 27 , 25 , 27 , 25 , ...
                  27 , 25 , 40 ];
%  ------------------------------------------------------------------------
mkdir('Seed_Data');
count = 0;
for i = 1:1:63
    Type = Lattice_Type_Matrix(i,:);
    Max_Variation = Variation_Matrix(i,:);
    if Max_Variation(1) ~= 0 && Max_Variation(2) == 0
        Mu_Vector = [0:Max_Variation(1)/4:Max_Variation(1)];
        Sigma_Vector = [0:Max_Variation(3)/4:Max_Variation(3)];
        for ii = 1:1:length(Mu_Vector)
            for jj = 1:1:length(Sigma_Vector)
                count = count + 1;
                fprintf('-- Generating %d of %d seed data ... \n', ...
                                                 count,sum(Sample_Number));
                Seed_Number = round((Max_Seed_Number-Min_Seed_Number) * ...
                                                     rand+Min_Seed_Number);
                Variation = [Mu_Vector(ii),0,Sigma_Vector(jj)];
                [ Seed ] = Seed_Generator ...
                               ( Box , Seed_Number , Type , Variation );
                cd('Seed_Data');
                File_Name = strcat('Seed_',num2str(count),'.txt');
                dlmwrite(File_Name,Seed);
                cd ..
            end
        end
    end
    if Max_Variation(1) ~= 0 && Max_Variation(2) ~= 0
        Mu_Vector = [0:Max_Variation(1)/2:Max_Variation(1)];
        k_Vector = [0:Max_Variation(2)/2:Max_Variation(2)];
        Sigma_Vector = [0:Max_Variation(3)/2:Max_Variation(3)];
        for ii = 1:1:length(Mu_Vector)
            for jj = 1:1:length(k_Vector)
                for kk = 1:1:length(Sigma_Vector)
                    count = count + 1;
                    fprintf('-- Generating %d of %d seed data ... \n', ...
                                                 count,sum(Sample_Number));
                    Seed_Number = round((Max_Seed_Number - ...
                                         Min_Seed_Number) * ...
                                         rand+Min_Seed_Number);
                    Variation = [Mu_Vector(ii) , ...
                                    k_Vector(jj) , ...
                                    Sigma_Vector(kk)];
                    [ Seed ] = Seed_Generator ...
                               ( Box , Seed_Number , Type , Variation );
                    cd('Seed_Data');
                    File_Name = strcat('Seed_',num2str(count),'.txt');
                    dlmwrite(File_Name,Seed);
                    cd ..
                end
            end
        end
    end
    if Max_Variation(1) == 0
        Sigma_Vector = [0:Max_Variation(3)/39:Max_Variation(3)];
        for ii = 1:1:length(Sigma_Vector)
            count = count + 1;
            fprintf('-- Generating %d of %d seed data ... \n', ...
                                             count,sum(Sample_Number));
            Seed_Number = round((Max_Seed_Number-Min_Seed_Number) * ...
                                                 rand+Min_Seed_Number);
            Variation = [0,0,Sigma_Vector(ii)];
            [ Seed ] = Seed_Generator ...
                           ( Box , Seed_Number , Type , Variation );
            cd('Seed_Data');
            File_Name = strcat('Seed_',num2str(count),'.txt');
            dlmwrite(File_Name,Seed);
            cd ..
        end
    end            
end
%  ------------------------------------------------------------------------
fprintf('-- All seed data have been generated! \n');
% =========================================================================