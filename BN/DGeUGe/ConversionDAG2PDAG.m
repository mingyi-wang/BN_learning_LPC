function cpdag=ConversionDAG2PDAG(dag)
% Conversion DAG -> CPDAG
% Adapted from code by Marco Grzegorczyk.
% Uses the functions topological_sort, parents, children
% by Kephin Murphy (from BNT)
% 
% INPUT: dag
% OUTPUT: cpdag 
% INVOCATION: cpdag=ConversionDAG2PDAG(dag) 

[N,M]=size(dag);
if N~=M
   error('DAG must be a square matrix');
end
cpdag = dag_to_cpdag(dag);

% cpdag ist eine Matrix, deren Einträge -1,0 oder 1 sind
% mit der Interpretation:
%  0 es gibt gar keine Kante
% +1 die Kante ist "compelled" (eindeutig festgelegte Orientierung)
% -1 die Kante ist "reversible" (keine eindeutige Orientierung)
X = (cpdag == 1);
Y = (cpdag == -1);

cpdag= X+Y+Y';

% -------------------------------------
function label = dag_to_cpdag(dag)
% -------------------------------------
%   Konvertierung DAG zu CPDAG
%       +1 falls Kante "compelled" 
%       -1 falss Kante "reversible"

N=length(dag);
[order xedge yedge] = order_edges(dag); % natuerliche Anordnung der Kanten
% Wichtig sind die geordneten Indices xedge und yedge

label = 2*dag; % Jede Kante wird mit "unknown" (=2) belegt

NbEdges = length(xedge) ; % Anzahl der Kanten


for Edge=1:NbEdges, % Fuer jede Kante:
    xlow=xedge(Edge); % Auswahl gemaess natuerlicher
    ylow=yedge(Edge); % Anordnung der Kanten
    
    if label(xlow,ylow)==2 % Falls die Kante unbekannt
        fin = 0; 
        wcompelled = find(label(:,xlow)==1); % Jede Kante, die nach xlow fuehrt und mit "compelled" belegt ist
        parenty = find(label(:,ylow)~=0); % Jede Kante, die nach ylow fuehrt (compelled oder unknown)
                                          % entspricht also allen Eltern von ylow

        for w = wcompelled' % Fuer alle Kanten nach xlow, die compelled sind
           
            if ~ismember(w,parenty) % Wenn w kein Elter von y
                label(parenty,ylow)=1; % Label xlow->ylow und alle w: w->ylow mit compelled (=1)
                fin = 1;  % Ende
            elseif fin == 0 % Sonst elseif-Anweisung:
                
                label(w,ylow)=1; % label w->y mit compelled (=1)-> fuehrt zu einem Fehler!!!!
               
            end
            
        end
        
        if fin == 0
            parentx = [xlow ; find(label(:,xlow)~=0)];
            % xlow und alle Kanten die nach xlow fuehren
            
            if ~isempty(mysetdiff(parenty,parentx)) % alle Eltern von y ausser x und den Eltern von x
                label(find(label(:,ylow)==2),ylow)=1; % label xlow->ylow und alle unknown Kante, die zu ylow fuehren mit compelled 
            else        
                label(xlow,ylow)=-1; % xlow->ylow "Reversible" 
                %label(ylow,xlow)=-1; % und umgekehrt
                label(find(label(:,ylow)==2),ylow)=-1; % Alle unknown-Kanten die nach ylow fuehren
                %label(ylow,find(label(:,ylow)==2))=-1; % werden "reversible" gesetzt und umgekehrt
            end
        end
    end
end

%%%========================================================================================
function [order, x, y] = order_edges(dag)
% ORDER_EDGES Algorithmus - ordne die Kanten in natuerlicher Weise 

N=length(dag);  % Anzahl der Knoten im Graphen
order = zeros(N,N); % Eine (NxN)-Matrix mit Nullen

% Topologische Sortierung der Knoten im Graphen:
node_order = topological_sort(dag);
[tmp oo] = sort(node_order);
% Es gilt: node_order = oo

% Insbesondere aber ist oo die Permutation, mit der topologisch sortiert wurde
% Es gilt: node_order(oo)=tmp UND tmp ist sortiert

dag=dag(node_order,node_order); % Die Adjacency-Matrix topologisch sortiert

% 1) flipud - Vertauschung der Zeilen (1te Zeile wird letzte Zeile usw.)
% 2) find - Bestimme die Koordinaten der "1"-Elemente 
[x y]=find(flipud(dag)==1);
% x und y sind die Koordinaten der Kanten, 
% nach topologischer Sortierung und Zeilenvertauschung

% Dies fuehrt zu einer Sortierung genau umgekehrt zur topologischen Sortierung
% Die letzte Zeile und die letzte Splate korrespondiert zu einem Knoten, der
% keine Eltern hat usw.

nb_edges=length(x); % Anzahl der Kanten

if nb_edges~=0  % Falls es Kanten gibt:
    % sub2ind liefert 1-dim Index der Eintraege (N+1-x,y) einer (N,N)-Matrix
    % Der Uebergang x -> N+1-x ist aufgrund der Zeilenvertauschung notwendig
  order(sub2ind([N N],N+1-x,y))=1:nb_edges ; % Die Kanten werden hier durchnummeriert
end

order=order(oo,oo); % Es wird zuruecksortiert, so dass wieder die Form von dag entsteht
% Aber in der order-Matrix sind die Kanten nun durchnummeriert

x=node_order(N+1-x);
y=node_order(y);
% x und y sind ihre Koordinaten entsprechend dieser Reihenfolge

% ---------------------------------------
function order = topological_sort(A)
% ---------------------------------------
% TOPOLOGICAL_SORT Return the nodes in topological order (parents before children).
% order = topological_sort(adj_mat)

n = length(A);
indeg = zeros(1,n);
zero_indeg = []; % a stack of nodes with no parents
for i=1:n
  indeg(i) = length(parents(A,i));
  if indeg(i)==0
    zero_indeg = [i zero_indeg];
  end
end

t=1;
order = zeros(1,n);
while ~isempty(zero_indeg)
  v = zero_indeg(1); % pop v
  zero_indeg = zero_indeg(2:end);
  order(t) = v;
  t = t + 1;
  cs = children(A, v);
  for j=1:length(cs)
    c = cs(j);
    indeg(c) = indeg(c) - 1;
    if indeg(c) == 0
      zero_indeg = [c zero_indeg]; % push c 
    end
  end
end

% ---------------------------------
function ps = parents(adj_mat, i)
% ---------------------------------
% PARENTS Return the list of parents of node i
% ps = parents(adj_mat, i)

ps = find(adj_mat(:,i))';


% -----------------------------------
function cs = children(adj_mat, i, t)
% -----------------------------------
% CHILDREN Return the indices of a node's children in sorted order
% c = children(adj_mat, i, t)
%
% t is an optional argument: if present, dag is assumed to be a 2-slice DBN

if nargin < 3 
  cs = find(adj_mat(i,:));
else
  if t==1
    cs = find(adj_mat(i,:));
  else
    ss = length(adj_mat)/2;
    j = i+ss;
    cs = find(adj_mat(j,:)) + (t-2)*ss;
  end
end

% -------------------------
function C = mysetdiff(A,B)
% -------------------------
% MYSETDIFF Set difference of two sets of positive integers (much faster than built-in setdiff)
% C = mysetdiff(A,B)
% C = A \ B = { things in A that are not in B }

if isempty(A)
  ma = 0;
else
  ma = max(A);
end

if isempty(B)
  mb = 0;
else
  mb = max(B);
end

if ma==0 
  C = [];
elseif mb==0
  C = A;
else % both non-empty
  %bits = sparse(1, max(ma,mb));
  bits = zeros(1, max(ma,mb));
  bits(A) = 1;
  bits(B) = 0;
  C = A(logical(bits(A)));
end
