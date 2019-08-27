function [NewBayes]=MCMCOCheckx1_x2(SCORES,NewBayes,neworder,order,x1,x2);

for pos=(x1+1):(x2-1)
    node=neworder(pos);
    nro_families=length(SCORES{node}.Parents);

    for family=1:nro_families
        flag_hasx1=0;
        flag_hasx2=0;
        flag_family_was_valid=1;
        flag_family_is_valid=1;

        nro_parents=length(SCORES{node}.Parents{family});
        if nro_parents>0;
            for parent_nro=1:nro_parents
                if SCORES{node}.Parents{family}(parent_nro)==neworder(x1)
                    flag_hasx1=1;
                end
                if SCORES{node}.Parents{family}(parent_nro)==neworder(x2)
                    flag_hasx2=1;
                end
            end%parents

            if flag_hasx1==1
                for parent_nro=1:nro_parents
                    position_parent=find(neworder==SCORES{node}.Parents{family}(parent_nro));
                    if position_parent>=pos
                        flag_family_is_valid=0;
                    end
                end%parents
                if flag_family_is_valid==1
                    %node 
                    %family
                    NewBayes(node,family)=NewBayes(node,family)+exp(SCORES{node}.LogBayesScores(family));
                end
            end%ifhasx1

            if flag_hasx2==1
                for parent_nro=1:nro_parents
                    position_parent=find(order==SCORES{node}.Parents{family}(parent_nro));
                    if position_parent>=pos
                        flag_family_was_valid=0;
                    end
                end%parents
                if flag_family_was_valid==1
                    NewBayes(node,family)=NewBayes(node,family)-exp(SCORES{node}.LogBayesScores(family));
                end
            end%ifhasx2

        end%ifpar>0

    end%families
end%pos