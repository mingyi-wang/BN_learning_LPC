function [NewBayes]=MCMCOCheckx1(SCORES,NewBayes,neworder,order,x1,x2)

node=neworder(x1);
nro_families=length(SCORES{node}.Parents);

for family=1:nro_families
    flag_par_betweenx1x2=0;
    flag_family_was_valid=1;
    nro_parents=length(SCORES{node}.Parents{family});
    if nro_parents>0;
        for parent_nro=1:nro_parents
            position_parent=find(neworder==SCORES{node}.Parents{family}(parent_nro));
            if position_parent>=x1 & position_parent<=x2
                flag_par_betweenx1x2=1;
                %family
            end
        end%parents
        if flag_par_betweenx1x2==1
            for parent_nro=1:nro_parents
                position_parent=find(order==SCORES{node}.Parents{family}(parent_nro));
                if position_parent>=x2
                    flag_family_was_valid=0;
                end
            end%parents
            if flag_family_was_valid==1
                NewBayes(node,family)=NewBayes(node,family)-exp(SCORES{node}.LogBayesScores(family));
                %family
            end%iffamiwas
        end%ifflagbet

    end%ifpar>0

end%families