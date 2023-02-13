%  Writting Abaqus Input File for Uniaxial Compression
%  Johns Hopkins University
%  Shengzhi Luan
%  02.08.2023
% =========================================================================
function  Abaqus_Input ( Material_Property )
%  ------------------------------------------------------------------------
%  Input:  Material_Property -- Material.
%  ------------------------------------------------------------------------
File = strcat('Lattice','.inp');
Introduction = 'Lattice Uniaxial Compression';
Company = 'Johns Hopkins University';
Author = 'Shengzhi Luan';
Data = '02.08.2023';
[ File ] = File_Introduction ( File , Introduction , ...
                               Company , Author , Data );
%  ------------------------------------------------------------------------
Heading = "Lattice Uniaxial Tension";
[ File ] = File_Heading ( File , Heading );
[ File ] = File_Double_Rule_Line ( File );
%  ------------------------------------------------------------------------
Part_Name = "Lattice";
[ File ] = File_Part_Creation ( File , Part_Name );
Node_Name = "Lattice";
[ File ] = File_Node ( File , Node_Name );
Element_Name = "Lattice";
Element_Type = "B22";
[ File ] = File_Element ( File , Element_Name , Element_Type );
%  ------------------------------------------------------------------------
Section_Type = "Beam";
Section_Element = "Lattice";
Material = ["Resin10K"];
Section_Shape = "Rect";
Thickness = Material_Property(1);
Section_Data = strcat(num2str(Thickness),', ',num2str(Thickness));
[ File ] = File_Section ( File , Section_Type , Section_Element , ...
                                 Material , Section_Shape , Section_Data );
[ File ] = File_Part_End ( File );
[ File ] = File_Double_Rule_Line ( File );
%  ------------------------------------------------------------------------
Assembly_Name = "Lattice";
[ File ] = File_Assembly_Creation ( File , Assembly_Name );
[ File ] = File_Single_Rule_Line ( File );
%  ------------------------------------------------------------------------
Part_Name = "Lattice";
Instance_Name = "Real_Lattice";
[ File ] = File_Instance ( File , Part_Name , Instance_Name );
%  ------------------------------------------------------------------------
load Bottom.txt;
load Top.txt;
Nset_Name = strcat(Instance_Name,'_Bottom');
[ File ] = File_Nset ( File , Nset_Name , Instance_Name , Bottom );
Nset_Name = strcat(Instance_Name,'_Top');
[ File ] = File_Nset ( File , Nset_Name , Instance_Name , Top );
[ File ] = File_Single_Rule_Line ( File );
%  ------------------------------------------------------------------------
[ File ] = File_Reference_Point ( File );
Node_1 = convertCharsToStrings(strcat(Instance_Name,'_Bottom'));
Node_2 = "Bottom_Ref";
Node_Vector = [Node_1;Node_2];
Coe_Vector = [1;-1];
DoF = [1,2,6];
[ File ] = File_Equation ( File , Node_Vector , Coe_Vector , DoF );
Node_1 = convertCharsToStrings(strcat(Instance_Name,'_Top'));
Node_2 = "Top_Ref";
Node_Vector = [Node_1;Node_2];
Coe_Vector = [1;-1];
DoF = [1,2,6];
[ File ] = File_Equation ( File , Node_Vector , Coe_Vector , DoF );
%  ------------------------------------------------------------------------
[ File ] = File_Assembly_End ( File );
[ File ] = File_Double_Rule_Line ( File );
%  ------------------------------------------------------------------------
Material_Name = "Resin10K";
General = ["Density";"1750"];
Mechanical = ["Elastic";"-";"10000E6, 0.30"];
[ File ] = File_Material ( File , Material_Name , General , Mechanical );
[ File ] = File_Double_Rule_Line ( File );
%  ------------------------------------------------------------------------
Step_Name = "Compression";
Step_Type = "Static";
Step_Data = ", , , ,";
Nlgeom_Switch = "No";
[ File ] = File_Step_Creation ( File , Step_Name , Step_Type , ...
                                               Step_Data , Nlgeom_Switch );
Boundary_Set = 'Bottom_Ref';
DoF = [1,2,6];
Boundary_Value = [0,0,0];
[ File ] = File_Boundary ( File , Boundary_Set , DoF , Boundary_Value );
Boundary_Set = 'Top_Ref';
DoF = [1,2,6];
Boundary_Value = [0,-1,0];
[ File ] = File_Boundary ( File , Boundary_Set , DoF , Boundary_Value );
[ File ] = File_Step_End ( File );
[ File ] = File_Double_Rule_Line ( File );
%  ------------------------------------------------------------------------
end
% =========================================================================



%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Assembly_Creation ( Input_File , ...
                                                     Assembly_Name )
%  Input:   Input_File -- Name of the input file.
%           Assembly_Name -- Name of the assembly.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Assembly, Name = ');
fprintf(File_ID,Assembly_Name);
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Assembly_End ( Input_File )
%  Input:   Input_File -- Name of the input file.
%           Assembly_Name -- Name of the assembly.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*End Assembly\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Boundary ( Input_File , Boundary_Set , ...
                                            DoF , Boundary_Value )
%  Input:   Input_File -- Name of the input file.
%           Boundary_Set -- Node(or element) set implied boundary.
%           DoF -- Degree of freedof of the boundary value.
%           Boundary_Value -- Displacement value of the boundary.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Boundary\n');
for i = 1:1:length(DoF)
    fprintf(File_ID,Boundary_Set);
    fprintf(File_ID,', ');
    fprintf(File_ID,'%d, %d, %e\n',DoF(i),DoF(i),Boundary_Value(i));
end
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Cload ( Input_File , Cload_Set , ...
                                         DoF , Cload_Value )
%  Input:   Input_File -- Name of the input file.
%           Boundary_Set -- Node(or element) set implied boundary.
%           DoF -- Degree of freedof of the boundary value.
%           Boundary_Value -- Displacement value of the boundary.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Cload\n');
for i = 1:1:length(DoF)
    fprintf(File_ID,Cload_Set);
    fprintf(File_ID,', ');
    fprintf(File_ID,'%d, %e\n',DoF(i),Cload_Value(i));
end
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Double_Rule_Line ( Input_File )
%  Input:   Input_File -- Name of the input file.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
number_of_line = 100;
number_of_comment = 2;
for i = 1:1:number_of_line
    Double_Rule_Line(i) = '=';
end
for i = 1:1:number_of_comment
    Double_Rule_Line(i) = '*';
end
for i = number_of_line-number_of_comment+1:1:number_of_line
    Double_Rule_Line(i) = '*';
end
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,Double_Rule_Line);
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Element ( Input_File , Element_Name , ...
                                           Element_Type )
%  Input:   Input_File -- Name of the input file.
%           Element_Name -- Name of the element.
%           Element_Type -- Type of the element.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Element, Type = ');
fprintf(File_ID,Element_Type);
fprintf(File_ID,', Input = T_');
fprintf(File_ID,Element_Name);
fprintf(File_ID,'.txt, Elset = ');
fprintf(File_ID,Element_Name);
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Elset ( Input_File , Elset_Name , ...
                                         Instance_Name , Element_ID )
%  Input:   Input_File -- Name of the input file.
%           Elset_Name -- Name of the elset.
%           Element_ID -- ID of the element.(Integer or integer vector)
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Elset, Elset = ');
fprintf(File_ID,Elset_Name);
fprintf(File_ID,', Instance = ');
fprintf(File_ID,Instance_Name);
fprintf(File_ID,'\n');
fprintf(File_ID,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',Element_ID);
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Equation ( Input_File , Node_Vector , ...
                                            Coe_Vector , DoF )
%  Input:   Input_File -- Name of the input file.
%           Node_Name -- Name of the node.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
for i = 1:1:length(DoF)
    fprintf(File_ID,'*Equation\n');
    fprintf(File_ID,num2str(length(Coe_Vector)));
    fprintf(File_ID,' \n');
    for j = 1:1:length(Coe_Vector)
        fprintf(File_ID,Node_Vector(j));
        fprintf(File_ID,', ');
        fprintf(File_ID,'%d',DoF(i));
        fprintf(File_ID,', ');
        fprintf(File_ID,'%f, ',Coe_Vector(j));
    end
    fprintf(File_ID,' \n');
end
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Heading ( Input_File , Heading )
%  Input:   Input_File -- Name of the input file.
%           Heading -- Heading of the input file.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Heading');
fprintf(File_ID,'\n');
fprintf(File_ID,Heading);
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Instance ( Input_File , Part_Name , ...
                                            Instance_Name )
%  Writing Instance into Abaqus Input File
%  Johns Hopkins University
%  By S.Z.Luan
%  03.30.2018
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Instance, Name = ');
fprintf(File_ID,Instance_Name);
fprintf(File_ID,', Part = ');
fprintf(File_ID,Part_Name);
fprintf(File_ID,'\n');
fprintf(File_ID,'*End Instance');
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Introduction ( Input_File , ...
                                   Introduction , Company , Author , Data )
%  Input:   Input_File -- Name of the input file.
%           Introduction -- Introduction of the input file.
%           Company -- Company of the input file.
%           Author -- Author of the input file.
%           Data -- Data of the input file.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
number_of_line = 100;
number_of_comment = 2;
for i = 1:1:number_of_line
    Star_Line(i) = '*';
    Introduction_Line(i) = ' ';
    Company_Line(i) = ' ';
    Author_Line(i) = ' ';
    Data_Line(i) = ' ';
end
if length(Introduction) <= number_of_line-2*2*number_of_comment
    for i = 1:1:2*number_of_comment
        Introduction_Line(i) = '*';
    end
    for i = number_of_line-2*number_of_comment+1:1:number_of_line
        Introduction_Line(i) = '*';
    end
    Introduction_Line(2*number_of_comment+2+1: ...
        2*number_of_comment+2+length(Introduction)) = Introduction;
end
if length(Introduction) > number_of_line-2*2*number_of_comment
    for i = 1:1:2*number_of_comment
        Introduction_Line(i) = '*';
    end
    for i = number_of_line-2*number_of_comment+1:1:number_of_line
        Introduction_Line(i) = '*';
    end
    Warning_Message = ...
        'Warning:The length of introduction is beyond length limitation.';
    Introduction_Line(2*number_of_comment+2+1: ...
        2*number_of_comment+2+length(Warning_Message)) = Warning_Message;
end
for i = 1:1:2*number_of_comment
    Company_Line(i) = '*';
end
for i = number_of_line-2*number_of_comment+1:1:number_of_line
    Company_Line(i) = '*';
end
Company_Line(2*number_of_comment+2+1: ...
    2*number_of_comment+2+length(Company)) = Company;
for i = 1:1:2*number_of_comment
    Author_Line(i) = '*';
end
for i = number_of_line-2*number_of_comment+1:1:number_of_line
    Author_Line(i) = '*';
end
Author_Line(2*number_of_comment+2+1: ...
    2*number_of_comment+2+length(Author)) = Author;
for i = 1:1:2*number_of_comment
    Data_Line(i) = '*';
end
for i = number_of_line-2*number_of_comment+1:1:number_of_line
    Data_Line(i) = '*';
end
Data_Line(2*number_of_comment+2+1: ...
    2*number_of_comment+2+length(Data)) = Data;
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,Star_Line);
fprintf(File_ID,'\n');
fprintf(File_ID,Star_Line);
fprintf(File_ID,'\n');
fprintf(File_ID,Introduction_Line);
fprintf(File_ID,'\n');
fprintf(File_ID,Company_Line);
fprintf(File_ID,'\n');
fprintf(File_ID,Author_Line);
fprintf(File_ID,'\n');
fprintf(File_ID,Data_Line);
fprintf(File_ID,'\n');
fprintf(File_ID,Star_Line);
fprintf(File_ID,'\n');
fprintf(File_ID,Star_Line);
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Material ( Input_File , ...
                                     Material_Name , General , Mechanical )
%  Input:   Input_File -- Name of the input file.
%           Material_Name -- Name of the material.
%           General -- General properties and data of the material.
%           Mechanical -- General properties and data of the material.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Material, Name = ');
fprintf(File_ID,Material_Name);
fprintf(File_ID,'\n');
fprintf(File_ID,'*');
fprintf(File_ID,General(1));
fprintf(File_ID,'\n');
fprintf(File_ID,General(2));
fprintf(File_ID,',\n');
fprintf(File_ID,'*');
fprintf(File_ID,Mechanical(1));
if Mechanical(2) ~= "-"
    fprintf(File_ID,', ');
    fprintf(File_ID,Mechanical(2));
end
fprintf(File_ID,'\n');
fprintf(File_ID,Mechanical(3));
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
%  ------------------------------------------------------------------------
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Node ( Input_File , Node_Name )
%  Input:   Input_File -- Name of the input file.
%           Node_Name -- Name of the node.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Node, Input = P_');
fprintf(File_ID,Node_Name);
fprintf(File_ID,'.txt, Nset = ');
fprintf(File_ID,Node_Name);
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Nset ( Input_File , Nset_Name , ...
                                        Instance_Name , Node_ID )
%  Input:   Input_File -- Name of the input file.
%           Nset_Name -- Name of the nset.
%           Instance_Name -- Name of the instance.
%           Node_ID -- ID of the node.(Integer or integer vector)
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Nset, Nset = ');
fprintf(File_ID,Nset_Name);
fprintf(File_ID,', Instance = ');
fprintf(File_ID,Instance_Name);
fprintf(File_ID,'\n');
fprintf(File_ID,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',Node_ID);
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Part_Creation ( Input_File , Part_Name )
%  Input:   Input_File -- Name of the input file.
%           Part_Name -- Name of the part.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Part, Name = ');
fprintf(File_ID,Part_Name);
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Part_End ( Input_File )
%  Input:   Input_File -- Name of the input file.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*End Part');
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Reference_Point ( Input_File )
%  Input:   Input_File -- Name of the input file.
%           Step_Name -- Name of the step.
%           Step_Type -- Type of the step.
%           Step_Date -- Data of the step.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Node\n');
fprintf(File_ID,'1, 0.000000, 0.000000, 0.000000\n');
fprintf(File_ID,'*Nset, Nset = Bottom_Ref, Internal\n');
fprintf(File_ID,'1\n');
fprintf(File_ID,'*Node\n');
fprintf(File_ID,'2, 0.000000, 0.000000, 0.000000\n');
fprintf(File_ID,'*Nset, Nset = Top_Ref, Internal\n');
fprintf(File_ID,'2\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Section ( Input_File , Section_Type , ...
                Section_Element , Material , Section_Shape , Section_Data )
%  Input:   Input_File -- Name of the input file.
%           Section -- Type of the section.
%           Section_Element -- Element assigned by the section.
%           Material -- Name of the material in the section.
%           Section_Shape -- Shape of the section geometry.
%           Section_Data -- Special date of the section geometry.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
if Section_Type == "Beam"
    File_ID = fopen(Input_File,'a+');
    fprintf(File_ID,'*');
    fprintf(File_ID,Section_Type);
    fprintf(File_ID,' Section, Elset = ');
    fprintf(File_ID,Section_Element);
    fprintf(File_ID,', Material = ');
    fprintf(File_ID,Material);
    fprintf(File_ID,', Section = ');
    fprintf(File_ID,Section_Shape);
    fprintf(File_ID,'\n');
    fprintf(File_ID,Section_Data);
    fprintf(File_ID,'\n');
    fprintf(File_ID,'0., 0., -1.');
    fprintf(File_ID,'\n');
    fclose(File_ID);
    Output_File = Input_File;
end
if Section_Type == "Beam General"
    File_ID = fopen(Input_File,'a+');
    fprintf(File_ID,'*');
    fprintf(File_ID,Section_Type);
    fprintf(File_ID,' Section, Elset = ');
    fprintf(File_ID,Section_Element);
    fprintf(File_ID,', Section = General\n');
    fprintf(File_ID,Section_Data);
    fprintf(File_ID,'\n');
    fprintf(File_ID,'0., 0., -1.');
    fprintf(File_ID,'\n');
    fprintf(File_ID,Material);
    fprintf(File_ID,'\n');
    fclose(File_ID);
    Output_File = Input_File;
end    
if Section_Type == "Solid"
    File_ID = fopen(Input_File,'a+');
    fprintf(File_ID,'*');
    fprintf(File_ID,Section_Type);
    fprintf(File_ID,' Section, Elset = ');
    fprintf(File_ID,Section_Element);
    fprintf(File_ID,', Material = ');
    fprintf(File_ID,Material);
    fprintf(File_ID,'\n');
    if Section_Data ~= "-"
        fprintf(File_ID,'%f',Section_Data);
        fprintf(File_ID,'\n');
    end
    if Section_Data == "-"
        fprintf(File_ID,',\n');
    end
    fclose(File_ID);
    Output_File = Input_File;
end
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Single_Rule_Line ( Input_File )
%  Input:   Input_File -- Name of the input file.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
number_of_line = 100;
number_of_comment = 2;
for i = 1:1:number_of_line
    Double_Rule_Line(i) = '-';
end
for i = 1:1:number_of_comment
    Double_Rule_Line(i) = '*';
end
for i = number_of_line-number_of_comment+1:1:number_of_line
    Double_Rule_Line(i) = '*';
end
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,Double_Rule_Line);
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Step_Creation ( Input_File , ...
                        Step_Name , Step_Type , Step_Data , Nlgeom_Switch )
%  Input:   Input_File -- Name of the input file.
%           Step_Name -- Name of the step.
%           Step_Type -- Type of the step.
%           Step_Date -- Data of the step.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Step, Name = ');
fprintf(File_ID,Step_Name);
fprintf(File_ID,', Nlgeom = ');
fprintf(File_ID,Nlgeom_Switch);
if Step_Type == "Perturbation"
    fprintf(File_ID,', Perturbation');
end
fprintf(File_ID,'\n');
fprintf(File_ID,'*');
fprintf(File_ID,Step_Type);
fprintf(File_ID,'\n');
fprintf(File_ID,Step_Data);
fprintf(File_ID,'\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================




%  ============================== Subroutine ==============================
function  [ Output_File ] = File_Step_End ( Input_File )
%  Input:   Input_File -- Name of the input file.
%  Output:  Output_File -- Name of the output file.
%  ------------------------------------------------------------------------
File_ID = fopen(Input_File,'a+');
fprintf(File_ID,'*Restart, Write , Frequency = 0\n');
fprintf(File_ID,'*Output, Field, Variable = Preselect\n');
fprintf(File_ID,'*Element Output, Directions = Yes\n');
fprintf(File_ID,'SF\n');
fprintf(File_ID,'*Output, History, Variable = Preselect\n');
fprintf(File_ID,'*Node Print, Nset = Bottom_Ref\n');
fprintf(File_ID,'RF2\n');
fprintf(File_ID,'*End Step\n');
fclose(File_ID);
Output_File = Input_File;
end
%  ========================================================================