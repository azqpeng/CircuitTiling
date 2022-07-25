% Alex Peng, Last Edited: 6/13/22
% NOTE: THe first indexed node is the ground and the source. matlab why
% must you index starting at one. Also second indexed node is the sink.

% TO DO:

%tCirc is an NxNx2 matrix, such that the first layer is the adjacency
%matrix with current weights, and the second layer is the "rotation system
%matrix" with element i,j corresponding to the order of j in i's
%neighbors.note that i,j not necc equal j,i
function driver()
%pCirc = [0 2/3 1.5 0; 2/3 0 1 1.5; 1.5 1 0 2/3; 0 1.5 2/3 0]
%cCirc = [0 2 3 0; 0 0 1 2; 0 0 0 3; 0 0 0 0]

rCirc = [0 1 2 0; 1 0 2 3; 1 2 0 3; 0 1 2 0]
pCirc = [0 1 1 0; 1 0 1 1; 1 1 0 1; 0 1 1 0]
cCirc = [0 1 1 0; 0 0 0 1; 0 0 0 1; 0 0 0 0]





circuitTiler(pCirc, cCirc, rCirc)
end

function circuitTiler(pCirc, cCirc, rCirc)
% Function takes the graph rotation  circuit matrix as input (assumes the 
% 1st and 2nd column/row represent the source and sink respectively) and 
% plots the matrix.
% Inputs: rCirc, the rotation adjacency matrix thing (ordering of nodes)
% cCirc, the current edge adjacency matrix thing
% Outputs: squareTiling, the square tiling (data type??? TO DO)

% Initialization
n = size(cCirc, 1);

% going through every node and putting it in? problem is it can't recursion

% I think BFS might be the best option? how the fuck do you make the thing
% TODO: Actual queue class
visited = false([n 1]);
queue = [];

queue(end+1) = 1;
visited(1) = true;

coord = [0 0];

while ~isempty(queue)    
    % dequeue
    row = queue(1);
    queue = queue(2:end);
    mord = n + 1;

    for i=1:n
        ordList = zeros(n);
        ord = rCirc(row, i);
        % Figure out what edges exist
        if ord ~= 0
            if visited(i) == false
                queue(end+1)= i;
                visited(i) = true;
                ordList(ord) = i;
            end
        end
        for col=1:n
            if col~=0
                w = cCirc(row, col);
                h = w*pCirc(row, col);
                %TO DO above line is inefficient
                rectangle('Position', [coord w h]);
                print = [coord w h];
                coord(1) = coord(1) + w;
                coord(2) = coord(2) + h;
            else
               break
            end
        end
    end
end
end

function showCircuit(pCirc)
% Takes the circuit's adjacency matrix, converts it to a digraph object and
% plots it. NOTE: The first element in matrix is the ground and the source.
% Inputs: adjacency matrix "p_circuit"
% Outputs: n/a
G = digraph(pCirc);
plot(G, 'EdgeLabel', G.Edges.Weight)
end

function cCirc = findCurrent(pCirc, voltpot)
% Uses MNA algorithm to find the currents through each edge. 
% Inputs: pCirc, n x n adjacency holding all the resistances
% voltpot, nx1 voltage difference between nodes and source.
% Output: cCirc, n x n adjacency holding all the currents 

% Adjusts voltpot to include source, initializes cCirc and visited
voltpot = circshift(voltpot, 1);
voltpot(1,1) = 0;
n = size(pCirc, 1);
cCirc = zeros(n);
visited = [];
% Does DFS but checks if edges been visited and if not calculates current.
mDFS(pCirc, voltpot, cCirc, 0, visited);
end

function mDFS(pCirc, voltpot, cCirc, start, visited)
% Runs DFS but checks if not visited; if not visited then adds current
% through edge to cCirc
visited(end+1) = start;
for i = 1:size(pCirc, 1)
    if pCirc(start, i)~=0
        if ~ismember(i, visited)
            cCirc(start, i) = (voltpot(i) - voltpot(start)) / pCirc(i, start);
            mDFS(pCirc, voltpot, cCirc, i, visited)
        end
    end
end
end

function voltpot = mna(pCirc)
% Implements MNA algorithm for a DC planar circuit and one voltage source
% Inputs: pCircuit, an adjacency matrix such that the weight of each point
% is the resistance of an edge
% Outputs: "x" is the output matrix and holds unknown voltages plus current
% through the one voltage source (n-1+1 x 1), 

% Initialize the matrices. 
V = 100;
n = size(pCirc,1);
% "A" holds conductace into various nodes
A = zeros([n n]);
% "z" holds n zeros -- representing current sources -- and  one value at 
% end to represent voltage force (n-1+1 x 1). ***CHANGE for eff
z = zeros([n 1]);
z(n, :) = V;

% Iterate through every element in pCircuit and put them into "A" (G)
for i = 1:n
    for j=1:n
        conductance = 1/pCirc(i, j);
        % Checks if it's element to itself, then skips iteration
        if i == j
            continue
        % Checks edges connected to ground node and adds appropriately
        elseif i==1 || j==1
            m = max([i j])-1;
            A(m, m) = A(m, m) + conductance;
        % Adds normal things
        else
            A(i-1, j-1) = - conductance;
            A(j-1, i-1) = - conductance;
            A(i-1, i-1) = A(i-1, i-1) + conductance;
        end
     end
end

% Adds in the voltage sink values
A(n, 1) = 1;
A(1, n) = 1;

% Output I guess
voltpot = A\z;
end

