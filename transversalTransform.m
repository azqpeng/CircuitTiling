function transversalTransform()
% left --> right
st_out = [11 11 11 1 1 2 2 3 4 4 5 6 7 8 9 10];
st_in = [1 2 3 4 5 5 6 9 7 8 8 9 12 10 12 12];
st_color = repmat([255, 0, 0], 16, 1);

% bottom --> top
ab_out = [13 13 3 3 9 9 2 6 6 1 5 8 10 4 7];
ab_in = [3 9 2 6 8 10 1 5 8 14 4 7 7 14 14];
ab_color = repmat([0, 0, 255], 15, 1);

% combine the representations and plot
all_out = horzcat(st_out, ab_out);
all_in = horzcat(st_in, ab_in);
inout_mat = [all_out; all_in]';
all_color = uint8([st_color; ab_color]);

transversal = digraph(all_out, all_in, table(all_color));

% THIS DOESN'T PRESERVE THE ORDER WE WANT. figure out how to fix this

plot(transversal, "EdgeColor", transversal.Edges.all_color);
end



function rect_loc = transversal(transversal, start, x_val)
% Function takes transversal structure and converts it into a rectangular 
% dissection
% Inputs: transversal -- graph with red and blue labels, start -- left
% node.

% Initialization
node_num = numnodes(transversal);
rect_loc = zeros(4, node_num - 4);
node_names = transversal.Nodes;
visited = containers.Map(node_names, false(size(node_names)));
queue = [];

queue(end+1) = start;
visited(start) = true;
coord = [0 0];
flag = 0;

while ~isempty(queue)
    % dequeue
    curr = queue(1);
    queue = queue(2:end);
    nbrs = successors(transversal, curr);
    first = 0;
    % iterates through blue neighbors
    for nbr = nbrs
       % BFS-ing
       if visited(nbr) == false
            visited(nbr) = true;
            queue(end+1) = nbr;
       end
       
       % checks colors
       if transversal.Edges.Color(findedge(transversal, curr, nbr)) == 255
            for check = predecessors(transversal, nbr)
                if ismember(check, nbrs)
                    flag = 1;
                    break
                else
                    flag = 0;
                end
            end
       end

       % If the nbr isn't the last one, yay! otherwise keep going
        if flag == 1
            flag = 0;
            continue
        else
            first = nbr;
        end
    end
    
    coord = rect_loc(curr) + x_val(curr);
    nbrs = nbrs(nbrs~=first);
    
    height = x_val(first);
    while ~isempty(nbrs)
        w = x_val(first);
        rect_loc(first) = [coord, w, w];
        % rect_loc update
        % 
        % 
        coord = coord + height;

        % updates nbrs
        nbrs = nbrs(nbrs~=first);
        % updates "first" and then exits the loop
        for new = successors(transversal, first)
            if ismember(new, nbrs)
                first = new;
            end
        end
    end
end
end

function constructDual(transversal, start)

coord = [0 0];
botPath = true;
node_names = transversal.Nodes;
node_num = numnodes(transversal);

visit = containers.Map(node_names, false(size(node_names)));

% initialize boundaries
transversal.Nodes.left = zeros(node_num);
transversal.Nodes.right = zeros(node_num);
transversal.Nodes.top = zeros(node_num);
transversal.Nodes.bot = zeros(node_num);

% does some placing
place(start);

% Rightmost boundaries of rightmost path
right_path = rightPath(transversal);
for v = right_path
    transversal.Nodes.right(v) = coord(1) + 1;
end
end

function place(v)

if v ~= start 
    visit(v) = true;
    transversal.Nodes.left(v) = x;
    transversal.Nodes.top(v) = y;
    if FirstPath
        y = y+1;
        transversal.Nodes.bot(v) = y;
        max_y = y;
    else
        transversal.Nodes.bot(v) = y;
        for u=successors(v)
            if transversal.Edges.Color(findedge(transversal, v, u)) == 255
                if visit(u) & transversal.Nodes.bot(v) ~=transversal.Nodes.top(v)
                    transversal.Nodes.right(u) = x;
                    if transversal.Nodes.bot(u) > transversal.Nodes.bot(v)
                        y = transversal.Nodes.bot(u) + max([transversal.Nodes.top(u) transversal.Nodes.top(v)]) / 2;
                        transversal.Nodes.bot(v) = y;
                    end
                end
            end
        end
    end
end
if findNode(transversal, "bottom") == 0
    disp("error: bottom doesn't exist and you messed up")
end
if ismember("bottom", successors(v))
    transversal.Nodes.bot(v) = max_y;
else
    FirstChild = true
    for i=1
    end
end
end







function right_path = rightPath(transversal)
% This function takes a transversal structure as input and finds it's
% rightmost path.
% Input: transversal (graph) - a transversal structure
% Output: right_path (array) - a 1 x n array in order of rightmost path
% such that n is the number of elements in path

% Finds unordered list of rightmost path
if findnode(transversal, "right")==0
    disp("No Node Found");
else
    right_nodes = predecessors(transversal, "right");
    right_path = zeros(size(right_nodes));
end

% Finds first node in the rightmost path
for node=right_nodes
    if ismember("start", predecessors(node))
        next = node;
        i=1;
        right_path(i) = next;
        break
    end
end

% Successively checks nodes in the rightmost path
while next~=0
    for pot=successors(transversal, next)
        if ismember(pot, right_nodes)
            next = pot;
            i = i+1;
            right_path(i) = next;
        end
    end
end
end






