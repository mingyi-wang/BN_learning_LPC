function flag_validParents=OrderMCMC_ValidParents(nodeOrder,node,parentConfig)
% Decide if a node-parent configuration is valid given an order
% Dirk Husmeier, December 2004
%
% INPUT
% node - the node
% parentConfig - a vetor of parents
% nodeOrder - an order of nodes
%
% OUTPUT
% flag_validParents
% 0 --> parents invalid
% 1 --> parents valid
%
% INVOCATION
% flag_validParents=OrderMCMC_validParents(nodeOrder,node,parentConfig)

% Find rank of node
nodeRank=find(nodeOrder==node);

% Initilialise flag to valid
flag_validParents=1;

% Switch flag if one of the parents appears
% behind the node in the given order
Nparents= length(parentConfig);
if Nparents>0
   for nParent=1:Nparents
       parent=parentConfig(nParent);
       parentRank=find(nodeOrder==parent);
       if parentRank>=nodeRank
          flag_validParents=0;
       end
   end
end
