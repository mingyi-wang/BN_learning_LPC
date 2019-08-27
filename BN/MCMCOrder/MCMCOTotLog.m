function [bayes]=MCMCOTotLog(SCORES,order)
%% Given the scores and the order is calculated the total score to this
%% configuration. The result is a matrix nodes x families that are
%% contributing to the total scores. to get the scores pass this matrix to
%% the sum_small

nro_nodes=length(SCORES);
for node=1:nro_nodes
    size_family(node)=length(SCORES{node}.Parents);
end
bayes=zeros(nro_nodes,max(size_family));

for node=1:nro_nodes
    nro_families=length(SCORES{node}.Parents);
    position_node=find(order==node);
    for family=1:nro_families
        nro_parents=length(SCORES{node}.Parents{family});
        flag_family_valid=1;
        if nro_parents>0
            for parent_nro=1:nro_parents
                position_parent=find(order==SCORES{node}.Parents{family}(parent_nro));
                if position_parent>=position_node
                   flag_family_valid=0;
                end
            end%parents
        end%if parents>0
        if flag_family_valid==1
           bayes(node,family)=exp(SCORES{node}.LogBayesScores(family));
        end
    end%family
  end%nodes
%size(bayes)