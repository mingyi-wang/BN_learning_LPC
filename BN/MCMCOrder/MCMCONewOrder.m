function [NewOrder,n1,n2]=MCMCONewOrder(order)
%% this funcion only change the position of two nodes, creating a new order
%% and informs in wich positions the nodes where exchanged


size_order=length(order);
nodes_to_change=randperm(size_order);
n1=(nodes_to_change(1));
n2=(nodes_to_change(2));
temp=order(n1);
order(n1)=order(n2);
order(n2)=temp;
NewOrder=order;