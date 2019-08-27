function SCORES=SBmcmc_MakeScoreStructureGaussian(varargin)
% Computes score structure in a preprocessing step
% Original code Dirk Husmeier, July 2004
% adapted Adriano V. Werhli, August 2005
% Like MarcoMI/LogLhoodContributions(data,alpha),
% this function computes all contributions to the 
% Bayes score from individual parent-child constellations 
% up to cardinality three.
% The function differs from MarcoMI/LogLhoodContributions
% in the following respect.
% - It computes the log Bayes scores rather than the
%   marginal log likelihoods, i.e., a prior can optionally
%   be included.
% - Pruning is applied, that is, only relevant configurations
%   are kept. A relevant configuration is defined as follows:
%   The score must be among the best Nbest scores and the
%   difference to the best score must not exceed a pre-specified 
%   threshold.
%
% INPUT, MANDATORY
% - DATA: Training data modified by interventions. (data_to_DATA.m)
% - data: Training data.
% Rows are genes, columns are experiments.
% T_0, T_M, v, alpha--> These are prior parameters
% - Nbest
% The number of best configurations to be stored.
% Default: Nbest=100
% - minDifference
% Configurations are stored until 
% logScoreBest-logScore > minDifference
% Default: minDifference=10
% - flag_prior
% flag_prior==1
% This prior assumes that all
% parent cardinalities are equally likely.
% Note that the prior over structures is uniform
% flag_prior==0
% Uniform prior over structure
% Default: flag_prior=0
% 
% OUTPUT
% SCORES{node}.LogBayesScore: 
% Vector of scores, sorted in descending order.
% Example
% SCORES{node}.LogBayesScore(2): runner-up score
% SCORES{node}.Parents: Corresponding parent configurations
% Example
% SCORES{node}.Parents{2}: runner-up parent configuration;
% can be of dimension up to max-fan = 3; e.g.
% [] (no parent), [2] (1 parent), [1 5] (2 parents), 
% [2 6 7] (3 parents)
%
% FUNCTION CALLS
% LogLhoodContributionChildParents
%
% INVOCATION
% SBmcmc_MakeScoreStructureGaussian(data,T_0,T_m,v,alpha,'Nbest',Nbest,'minDifference',minDifference,'flag_prior',flag_prior);


% -----------------------------------
% Input
% -----------------------------------
DATA=varargin{1};
data=varargin{2};
T_0=varargin{3};
T_m=varargin{4};
v=varargin{5};
alpha=varargin{6};
[Nnodes,Ndata]=size(data);
%alpha=1;
Nbest=100;
minDifference=10;
flag_prior=0;
nargs = length(varargin);
for i=7:2:nargs
   switch varargin{i},
      case 'Nbest'
           Nbest=varargin{i+1};
      case 'minDifference'
           minDifference=varargin{i+1};
      case 'flag_prior'
           flag_prior=varargin{i+1};
      otherwise
            error('Wrong argument, possibly a mistyped keyword'); 
   end
end

% Number of parent configurations with different
% cardinalities
NParentsCard(1)= Nnodes-1;
NParentsCard(2)= (Nnodes-1)*(Nnodes-2)/2;
NParentsCard(3)= (Nnodes-1)*(Nnodes-2)*(Nnodes-3)/6;
NumberOfPossibleParentConfigurations=...
1 + ... % no parents
NParentsCard(1) + ... % 1 parent
NParentsCard(2) + ... % 2 parents
NParentsCard(3); % 3 parents
for i=1:3
    logNParentsCard(i)=log(NParentsCard(i));
end
if Nbest>NumberOfPossibleParentConfigurations
    error('The number of best configurations should be less than the number of possible parent configurations.\nWith this data set the number of possile configurations is %d',NumberOfPossibleParentConfigurations);
end
% =============================
%  Loop over nodes
% =============================
for node=1:Nnodes
LogBayesScores= zeros(1,NumberOfPossibleParentConfigurations);
% Vector of all the Bayes scores for a given node
clear ParentConfigs
% Structure with all the parent configurations
% for a given node.


node % Indicate progress on screen
nParent=0;

% -----------------------------------
% No parents
% -----------------------------------
nParent=nParent+1;
LogBayesScores(nParent)=...%%%%%%%%%%%%%%%%%%%%CHANGE%%%%%%%%%%%%%%%%%%%%%
Gauss_Score_complete_local_Intervention(node,[]',DATA,data,T_0,T_m,v,alpha);
%LogLhoodContributionChildParent(data,node,[],alpha);
ParentConfigs{nParent}=[];

% -----------------------------------
% 1 parent
% -----------------------------------
for parent1=1:Nnodes
    if parent1~=node
       nParent=nParent+1;
       LogBayesScores(nParent)=...%%%%%%%%%%%%%%CHANGE%%%%%%%%%%%%%%%%%%%%
       Gauss_Score_complete_local_Intervention(node,parent1,DATA,data,T_0,T_m,v,alpha);
       %LogLhoodContributionChildParent(data,node,parent1,alpha);
       ParentConfigs{nParent}=[parent1];
       if flag_prior
          LogBayesScores(nParent)=...
          LogBayesScores(nParent)-logNParentsCard(1);
       end
    end
end

% -----------------------------------
% 2 parents
% -----------------------------------
% Loop through all parent configurations for which
% each parent is different from the node. 
% Use the ordering parent1<parent2 to avoid double counts.
for parent1=1:Nnodes-1
    for parent2=parent1+1:Nnodes
        if parent1~=node & parent2~=node
           nParent=nParent+1;
           LogBayesScores(nParent)=...%%%%%%%%%%%%%%%%%%CHANGE%%%%%%%%%%%%%%%%%%%%%
           Gauss_Score_complete_local_Intervention(node,[parent1,parent2]',DATA,data,T_0,T_m,v,alpha);
           %LogLhoodContributionChildParent(data,node,[parent1,parent2],alpha);
           ParentConfigs{nParent}=[parent1,parent2];
           if flag_prior
              LogBayesScores(nParent)=...
              LogBayesScores(nParent)-logNParentsCard(2);
           end
        end
    end
end
% -----------------------------------
% 3 parents
% -----------------------------------
% Loop through all parent configurations for which
% each parent is different from the node. 
% Use the ordering parent1<parent2<parent3 to avoid double counts.
for parent1=1:Nnodes-2
    for parent2=parent1+1:Nnodes-1
        for parent3=parent2+1:Nnodes
            %parent1,parent2,parent3
            if parent1~=node & parent2~=node & parent3~=node
               nParent=nParent+1;
              % parent3
               LogBayesScores(nParent)=...%%%%%%%%%%%%%%%%%%%%%CHANGE%%%%%%%%%%%%%%%%%%%%
               Gauss_Score_complete_local_Intervention(node,[parent1,parent2,parent3]',DATA,data,T_0,T_m,v,alpha);
               %LogLhoodContributionChildParent(data,node,[parent1,parent2,parent3],alpha);
               ParentConfigs{nParent}=[parent1,parent2,parent3];
               if flag_prior
                  LogBayesScores(nParent)=...
                  LogBayesScores(nParent)-logNParentsCard(3);
               end
            end
        end
    end
end

% -----------------------------------
% Sorting and pruning
% -----------------------------------
[bestLogBayesScores,topIndices]= sort(-LogBayesScores);
bestLogBayesScores= -bestLogBayesScores;

% First step: Get the highest score
n=1;
SCORES{node}.LogBayesScores= bestLogBayesScores(n);
SCORES{node}.Parents{n}= ParentConfigs{topIndices(n)};

% Second step: Extend the list until the difference
% between the best and the worst configuration
% exceeds a pre-specified margin or until a total
% threshold number is exceeded
while (bestLogBayesScores(1)-bestLogBayesScores(n)) < minDifference & n<Nbest
  n=n+1;
  SCORES{node}.LogBayesScores(n)= bestLogBayesScores(n);
  SCORES{node}.Parents{n}= ParentConfigs{topIndices(n)};
end

end % -> node


% -----------------------------------
% Final output
% -----------------------------------
% SCORES;
