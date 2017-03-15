%% LQR-RRT* - Simple Pendulum
function pend_rrt( time_interval, N )

% LQR-RRT* is a fast motion planning algorithm which guarantees to find an 
% optimal solution asymptotically. Here the approach is used to diverge
% optimal motion policies for a simple pendulum in its phase plot. 
%
% Author: Mahan Fathi 
%
% Tree data is mainly stored in two matrices, one representing node
% properites and one edges. You don't need graph libraries on matlab to run
% this code.
%
% Notes:
%   1. Feel free to apply LQR-RRT* by using this as a template. 
%   2. Refinements such as Tree Prunning (Branch-and-Bound) or goal directed
%      propogation algorithms are yet to be implemented here and I'm looking 
%      forward for you to commit them.
%
% Inputs: 
%   1. time_interval: Determines how big you want your steps to be in time.
%   2. N            : Number of iterations the main loop counts.
%
% Output: 
%      Phase plot updates peridically in figure 25.
%

global GNodes GEdges model nun;           % nun --> Current Number of Nodes
GNodes  = inf * ones( N   , 2+1 );        % Each Row Contains States and Costz
GEdges  = inf * ones( N+1 , 2+1 );        % Each Row Contains EndNodes and Action

GNodes( 1, : )  = zeros( 1, 3 );
nun             = 1;

model.h                 = time_interval;
model.N                 = N;              % Number of Nodes

model.phy.m   = 1;
model.phy.l   = 1;
model.phy.b   = 0.1;
model.phy.g   = 9.81;

model.cost.Q   =  1*eye(2);
model.cost.R   = 50*eye(1);

model.bound.x1 = [ -pi pi ];
model.bound.x2 = [ -2.5*pi 2.5*pi ];
model.bound.u  = [ -5 5 ];


for i = 1 : model.N

SteerAction     = inf;
while abs(SteerAction) > 5
    x_rand = [ model.bound.x1(1)+(model.bound.x1(2)-model.bound.x1(1))*rand, ...
               model.bound.x2(1)+(model.bound.x2(2)-model.bound.x2(1))*rand ];
    % Setting Up Local LQR Cost and Policy on x_rand
    [ A_rand, B_rand ] = linmats( x_rand, 0 );
    [ K_rand, S_rand ] = lqr( A_rand, B_rand, model.cost.Q, model.cost.R );
    % Find Nearest Point in Phase Space
    [ id_nearest, x_nearest ]   = LQRNearest( x_rand, S_rand );
    % Construct SteerAction and Check If It's in Limits
    SteerAction        = -K_rand*( x_nearest - x_rand )';
end
% Steer from Nearest to New
x_new   = LQRSteer( x_nearest, x_rand, K_rand );
% Setting Up Local LQR Cost and Policy on x_new
[ A_new, B_new ] = linmats( x_new, 0 );
[ K_new, S_new ] = lqr( A_new, B_new, model.cost.Q, model.cost.R );
% Find the new neighborhood
X_near_ids = LQRNear( x_new, S_new, id_nearest );
% Choose Parent Candidate
[ id_parent, x_parent, minCost ] = ChooseParent( X_near_ids, x_new, S_new );
% If it's collision free add node x_new and the edge between parent
% canddidate and the new node in this case -- simple pendulum -- all new
% edges are collision free.
% New Node
to_parentCost     = GNodes( id_parent, 3 );
to_newCost        = to_parentCost + (x_parent-x_new)*S_new*(x_parent-x_new)';
nun               = nun + 1;
GNodes( nun, : )  = [ x_new to_newCost ];
% New Edge
action            = -K_new*( x_parent - x_new )';
GEdges( nun-1, : )= [ id_parent nun action' ];
% Rewire
Rewire( X_near_ids, x_new, i+1 );
% Draw the Graph

% Display Info
if mod(i,50)    ==  1
    fprintf('Iteration     EdgeCount    NodeCount\n')
    draw_graph( GNodes, GEdges )
end
fprintf( '  %i              %i            %i\n',i,nun-1,nun )

% Prune Costly Nodes
% if mod(i,40)    ==  1
%     fprintf('\n                 PRUNE\n\n')
%     Prune();
% end

% Check If Nodes Exists in Goal Neighborhood
% Check( G );

end

end



%===================================
%            Linearization
%===================================
function [ A, B ] = linmats( x0, u0 )

global model;
m   = model.phy.m;
l   = model.phy.l;
b   = model.phy.b;
g   = model.phy.g;

A   = [ 0, 1; -g*cos(x0(1))/l, -b/m/l^2 ];
B   = [ 0; 1/m/l^2 ];

end

%===================================
%            LQR Nearest
%===================================
function [ id_nearest, x_nearest ] = LQRNearest( x_rand, S_rand )

global GNodes nun;

nearest_cost    = inf;

    for i = 1 : nun
       x    = GNodes( i, 1:2 );
       cost = (x-x_rand)*S_rand*(x-x_rand)';
       if cost < nearest_cost
           nearest_cost     = cost;
           x_nearest        = x;
           id_nearest       = i;
       end
    end

end

%===================================
%            LQR Steer
%===================================
function x_new = LQRSteer( x_nearest, x_rand, K_rand )
% This function contains system's explicit dynamics
global model;
m   = model.phy.m;
l   = model.phy.l;
b   = model.phy.b;
g   = model.phy.g;
h   = model.h;

u       = -K_rand*( x_nearest - x_rand )';
x_dot   = [ x_nearest(2), ( u - b*x_nearest(2) - m*g*l*sin(x_nearest(1)) ) ];

x_new   = x_nearest + x_dot*h;

end

%===================================
%            LQR Near
%===================================
function X_near_ids = LQRNear( x_new, S_new, id_nearest )

global GNodes nun;

% Define Neighborhood Radius
gamma   = 1;  d   = 2;
ner     = gamma*( log(nun+1)/nun )^(1/d);
% Allocate and Assign to Output
X_near_ids  = id_nearest;

  for i = 1 : nun
    x    = GNodes( i, 1:2 );
    cost = (x-x_new)*S_new*(x-x_new)';
    if cost < ner
        X_near_ids(end+1)    = i;
    end
  end

X_near_ids  = unique(X_near_ids);

end

%===================================
%            Choose Parent
%===================================
function [ id_parent, x_parent, minCost ] = ChooseParent( X_near_ids, x_new, S_new )

global GNodes;

minCost     = inf;

    for id_candidate = X_near_ids
        x_candidate         = GNodes( id_candidate, 1:2 );
        to_canCost          = GNodes( id_candidate,  3  );
        from_canCost        = (x_candidate-x_new)*S_new*(x_candidate-x_new)';
        if to_canCost+from_canCost < minCost
            minCost     = to_canCost + from_canCost;
            id_parent   = id_candidate;
            x_parent    = x_candidate;
        end
    end

end

%===================================
%               Rewire
%===================================
function Rewire( X_near_ids, x_new, id_newnode )

global model GNodes GEdges;   % nun-1 --> Number of Edges

% Cost-to-go for x_new
to_newCost  = GNodes( id_newnode, 3 );

    for id_candidate = X_near_ids
        x_candidate     = GNodes( id_candidate, 1:2 );
        to_canCost      = GNodes( id_candidate,  3  );
        % Find Local System at x_near
        [ A_near, B_near ]  = linmats( x_candidate, 0 );
        [ K_near, S_near ]  = lqr( A_near, B_near, model.cost.Q, model.cost.R );
        ntocCost            = ( x_new-x_candidate )*S_near*( x_new-x_candidate )';
            if to_newCost + ntocCost < to_canCost
                % New Cost to Candidate Nodes Should be Modified in GNodes
                rnewcost  = to_newCost + ntocCost;
                GNodes( id_candidate, 3 ) = rnewcost;
                action      = -K_near*( x_new - x_candidate )';
                edge_index  =  GEdges(:,2) == id_candidate ;
                GEdges( edge_index, : )   = [ id_newnode id_candidate action ];
            end
    end

end

%===================================
%           Branch-and-Bound
%===================================
function Prune()

global G;

j = 1;
while j <= numnodes(G)
    if G.Nodes.Costz(j) > 5000
        G   = rmnode( G, j );
        j   = j - 1;
    end
    j = j + 1;
end


end

%===================================
%         Check Completeness
%===================================
function Check( G )






end


%===================================
%             Draw Graph
%===================================
function draw_graph( GNodes, GEdges )

global nun

persistent gFig
if (isempty(gFig))
    gFig = figure(25);
end

figure(gFig); cla; hold on;

for ids = GEdges( 1:nun-1, 1:2 )'

    plot( [GNodes(ids(1),1) GNodes(ids(2),1)], ...
          [GNodes(ids(1),2) GNodes(ids(2),2)] )

end

axis([ -pi pi -2.5*pi 2.5*pi ]);
drawnow

end
