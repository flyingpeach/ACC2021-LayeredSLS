function [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = setup_topology(gridSize, connectThresh)
% Set up square power grid without plant dynamics
%
% gridSize     : the square grid will be gridSize * gridSize
% connectThresh: connect neighbors in the grid with probability connectThresh
%                note that if a node is not connected to any neighbors we
%                will force connect it to one random neighbor
%
% adjMtx       : adjacency matrix of grid
% nodeCoords   : (x,y) coordinates of each node in 2D plane (for visualization purposes)
% susceptMtx   : matrix of susceptance (imaginary admittance) values between nodes
%

numNodes   = gridSize * gridSize;
nodeCoords = zeros(numNodes, 2);
for i=1:gridSize
    theseRows = (i-1)*gridSize+1:i*gridSize;
    nodeCoords(theseRows, 1) = i*ones(gridSize, 1);
    nodeCoords(theseRows, 2) = 1:gridSize;
end 

% set up adjacency matrix (randomly)
adjMtx        = zeros(numNodes, numNodes);

if numNodes > 1
    for thisNode = 1:numNodes
        leftNode   = thisNode - gridSize;
        rightNode  = thisNode + gridSize;
        topNode    = thisNode + 1;
        bottomNode = thisNode - 1;

        viableNodes = [leftNode, rightNode, topNode, bottomNode];

        if mod(thisNode, gridSize) == 1 % on the bottom
            viableNodes = viableNodes(viableNodes ~= bottomNode);
        elseif mod(thisNode, gridSize) == 0 % on the top
            viableNodes = viableNodes(viableNodes ~= topNode);
        end

        if leftNode < 1 % on the left
            viableNodes = viableNodes(viableNodes ~= leftNode);
        elseif rightNode > numNodes % on the right
            viableNodes = viableNodes(viableNodes ~= rightNode);
        end

        connected = false;
        % connect probabilistically
        for node=viableNodes
            coin = rand();
            if coin > connectThresh
                connected = true;
                adjMtx(thisNode, node) = 1; % symmetrical
                adjMtx(node, thisNode) = 1;            
            end
        end
    
        if connected == false % connect to at least one neighbour
            node = viableNodes(randi(length(viableNodes)));
            adjMtx(thisNode, node) = 1;
            adjMtx(node, thisNode) = 1;
        end 
    end
end

susceptMtx = zeros(numNodes, numNodes);

% make symmetrical
for i=1:numNodes
    for j=1:i
        if (adjMtx(i,j))
            susceptMtx(i,j) = rand()*0.5 + 0.5;
            susceptMtx(j,i) = susceptMtx(i,j);
        end
    end
end

% generate other physical parameters
inertiasInv = rand(numNodes, 1) * 10; % 1 / inertia
dampings    = zeros(numNodes, 1); % currently undamped

end