%  Display Tessellation from Tessellation Dataset
%  Johns Hopkins University
%  Shengzhi Luan
%  02.08.2023
% =========================================================================
Index = randi(1646);    %Input: Index of microstructure between 1 and 1646.
%  ------------------------------------------------------------------------
cd('Tessellation_Dataset');
Node = load(strcat('Node_',num2str(Index),'.txt'));
Connection = load(strcat('Connection_',num2str(Index),'.txt'));
cd ..
%  ------------------------------------------------------------------------
Frame = [0.15,0.85,0.15,0.85];
figure
grid off
hold on
for i = 1:1:size(Connection,1)
    plot(Node(Connection(i,2:3),2),Node(Connection(i,2:3),3), ...
         'Blue','Linewidth',1.5);
end
plot(Frame([1,2,2,1,1,2]),Frame([3,3,4,4,3,3]),'Black','Linewidth',2.0);
Title = strcat('Display of Tessellation: #',num2str(Index));
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