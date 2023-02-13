%  Quadratic Beam Element Discretization
%  Johns Hopkins University
%  Shengzhi Luan
%  02.08.2023
% =========================================================================
function  [ p_Output , t_Output ] = Beam_Discretization_Quadrature ...
                                          ( p_Input , t_Input , n_Element )
%  ------------------------------------------------------------------------
%  Input:   p_Input -- Nodal matrix ( 1st row is x-coordinate,
%                                     2nd row is y-coordinate,
%                                     3rd row is z-coordinate. ).
%           t_Input -- Connection matrix ( 1st row is start node ID,
%                                          2nd row is end node ID. ).
%           n_Element -- Number of element discretized in connection.
%  Output:  p_Output -- Nodal matrix ( 1st row is x-coordinate,
%                                      2nd row is y-coordinate,
%                                      3rd row is z-coordinate. ).
%           t_Output -- Element matrix ( 1st row is start node ID,
%                                        2nd row is middle node ID,
%                                        3rd row is end node ID. ).
%  ------------------------------------------------------------------------
p = p_Input;
t = t_Input;
n = n_Element;
Frame_Node_Number = size(p,2);
%  ------------------------------------------------------------------------
for i = 1:1:size(t,2)
    Node_A = p(:,t(1,i));
    Node_B = p(:,t(2,i));
    p_Discrete(:,(i-1)*(2*n-1)+1:(i-1)*(2*n-1)+n-1) = ...
               Node_A+[1/n:1/n:1-1/n].*(Node_B-Node_A);
    p_Discrete(:,(i-1)*(2*n-1)+n:(i-1)*(2*n-1)+2*n-1) = ...
               Node_A+[1/n/2:1/n:1-1/n/2].*(Node_B-Node_A);
    t_Discrete(:,(i-1)*n+1:(i-1)*n+n) = ...
        [0,(i-1)*(2*n-1)+1:(i-1)*(2*n-1)+n-1;
        (i-1)*(2*n-1)+n:(i-1)*(2*n-1)+2*n-1
        (i-1)*(2*n-1)+1:(i-1)*(2*n-1)+n-1,0] + Frame_Node_Number;
    t_Discrete(1,(i-1)*n+1) = t(1,i);
    t_Discrete(3,(i-1)*n+n) = t(2,i);
end
%  ------------------------------------------------------------------------
p = [p,p_Discrete];
t = [t_Discrete];
%  ------------------------------------------------------------------------
p_Output = p;
t_Output = t;
%  ------------------------------------------------------------------------
end
% =========================================================================