%  Display Microstructure from Tessellation Dataset
%  Johns Hopkins University
%  Shengzhi Luan
%  02.08.2023
% =========================================================================
Index = randi(1646);    %Input: Index of microstructure between 1 and 1646.
Density = randi(25)/100;       %Input: Relative density between 1% and 25%.
%  ------------------------------------------------------------------------
cd('Tessellation_Dataset');
Node = load(strcat('Node_',num2str(Index),'.txt'));
Connection = load(strcat('Connection_',num2str(Index),'.txt'));
cd ..
%  ------------------------------------------------------------------------
Sum_Length = sum(vecnorm(Node(Connection(:,2),[2,3])' - ...
                         Node(Connection(:,3),[2,3])'));
%  ------------------------------------------------------------------------
Frame = [0.15,0.85;0.15,0.85];
Frame_Area = prod(diff(Frame'));
Rho = 0;
count = 0;
Thickness = Frame_Area*Density/Sum_Length;
while abs(Rho-Density) >= 1E-5*Density && count <= 100
    for i = 1:1:size(Connection,1)
        Line = Polygon_Line ( Node(Connection(i,2),[2,3])' , ...
                              Node(Connection(i,3),[2,3])' , Thickness );
        if i == 1
            Metamaterial = Line;
        else
            Metamaterial = union(Metamaterial,Line);
        end
    end
    Rho = Metamaterial.area/Frame_Area;
    Thickness = Thickness*Density/Rho;
    count = count + 1;
    fprintf('Iteration #%d: relative density is %f ...\n',count,Rho);
end
%  ------------------------------------------------------------------------
figure
plot(Metamaterial,'FaceColor',[255,0,24]/255,'FaceAlpha',0.85, ...
                  'EdgeColor','Black','Linewidth',0.25,'EdgeAlpha',0.85);
grid off
hold on
plot(Frame(1,[1,2,2,1,1,2]),Frame(2,[1,1,2,2,1,1]),'Black','Linewidth',2);
Title = strcat('Display of Microstructure: #',num2str(Index), ...
               ' (',num2str(Density*100),'%)');
text(0.5,0.93,Title,'Color','Black', ...
                     'FontName','Times','FontSize',14, ...
                     'HorizontalAlignment','Center', ...
                     'VerticalAlignment','Middle');
axis([0.1,0.9,0.1,0.9]);
axis equal
set(gca,'XColor', 'none','YColor','none');
%  ------------------------------------------------------------------------
clear all
% =========================================================================




% =========================================================================
function  [ Line ] = Polygon_Line ( Node_A , Node_B , Thickness )
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
% =========================================================================