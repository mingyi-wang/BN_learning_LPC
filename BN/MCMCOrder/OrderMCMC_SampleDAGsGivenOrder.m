function [DAGs,logDAGscores]=OrderMCMC_SampleDAGsGivenOrder(ScoresStructure,nodeOrders,NsamplesPerOrder)
% Samples DAGs for a given set of orders
%
% INPUT
% ScoresStructure:  Structure of pre-computed log scores
%                   for node-parent configurations
% nodeOrders:       Matrix of node orders. 
%                   Each row is a different order.
% NsamplesPerOrder: Number of sampled DAGs per order.
%
% OUTPUT
% DAGs: Trajectory of DAGs
% logDAGscores:  Trajectory of log DAG scores
%
% FUNCTION CALL
% --> OrderMCMC_ScoresUnderOrderConstraint
% 
% INVOCATION
% [DAGs,logDAGscores]=OrderMCMC_SampleDAGsGivenOrder(ScoresStructure,nodeOrder)
tt=cputime;
% Get number of nodes
Nnodes=max(nodeOrders(1,:));
% Get number of sampled orders
Norders= size(nodeOrders,1);

ndag=0;
% --> Loop over all orders
for nOrder=1:Norders
    nOrder
    % Get node order
    nodeOrder=nodeOrders(nOrder,:);
    ScoresStructureUnderOrderConstraint=OrderMCMC_ScoresUnderOrderConstraint(ScoresStructure,nodeOrder);
    % --> Loop over repeated samples given an order
    for nSample=1:NsamplesPerOrder
      ndag=ndag+1;
      DAG= zeros(Nnodes); 
      totalNetScore=0;
      % --> Loop over all nodes: sample parent configurations
      for node=1:Nnodes
          node
          if node==80
             node 
          end 
	[parentConfig,localLogScore]= ProposalMove(node,ScoresStructureUnderOrderConstraint);
        totalNetScore=  totalNetScore+localLogScore;
        NumberOfParents= length(parentConfig);
        % --> Loop over all parents: create edges
        if NumberOfParents>0
           for nParent=1:NumberOfParents
	       parent=parentConfig(nParent);
               DAG(parent,node)=1;
           end
        end 
      end
      DAGs{ndag}=DAG;
      logDAGscores(ndag)= totalNetScore;
    end
end

ttotal=cputime-tt
% ----------------------------------------------------------------
    function [parentConfig,localLogScore]= ProposalMove(node,SCORES)
% ----------------------------------------------------------------
% Proposal move: New parent configuration 
% for a specified node.
% Importance sampling step: Propose new parents
% Ex.: CumProb=[0.1 0.4 0.6 0.9 1.0], rand=0.3
% proposed config = 2 (spanning the range 0.1->0.4). 
% However, for numerical stability the scores in SCORES
% are monotonically decreasing rather than increasing:
% CumProb=[1.0 0.9 0.6 0.4 0.1], rand=0.3
% proposed config = 4 (spanning the range 0.1->0.4). 
% This is achieved by the command:
% proposedID= sum(CumProb>=rand)
chain=1;
parentConfigID= sum(SCORES{node}.CumProb(chain,:)>=rand);
parentConfig= SCORES{node}.Parents{parentConfigID};
localLogScore=SCORES{node}.LogBayesScores(parentConfigID);
