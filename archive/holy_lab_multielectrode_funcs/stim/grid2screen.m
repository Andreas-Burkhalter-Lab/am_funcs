function screen_coordinates = grid2screen(stim_centers,gridColumnIndex,gridRowIndex)
%GRID2SCREEN: convert stim grid index into pixel coordinates on screen.
% Input:
%%% 1. stim_centers = cell array where stim_centers{row,column}(x,y) lists
%%%     the x,y coordinates of the stim center in the grid row index 'row'
%%%     and grid column index 'column' (L->R,top->bottom).
%%% 2. gridColumnIndex = grid column index (can be fractional) for which to
%%%     get the corresponding screen pixel x coordinate
%%% 3. gridRowIndex = grid row index (can be fractional) for which to
%%%     get the corresponding screen pixel y coordinate
% Output:
%%% screen_coordinates(x,y) = x,y coordinates on the screen in pixels
%%%     corresponding to the inputs gridColumnIndex and gridRowIndex, using
%%%     the coordinate conventions of stim_centers

% get pixel coordinates of the grid vertices to the left, right, above, and below
% the grid location of interest(usually will be between four grid vertices)
x_gridVertexToLeft = stim_centers{1,floor(gridColumnIndex)}(1);
x_gridVertexToRight = stim_centers{1,ceil(gridColumnIndex)}(1);
y_gridVertexAbove = stim_centers{floor(gridRowIndex),1}(2);
y_gridVertexBelow = stim_centers{ceil(gridRowIndex),1}(2);

% get displacement of screen_coordinates from the lower-index vertex 
% (to right or above) as a fraction of the grid spacing
x_extraFractionOfGridSpace = rem(gridColumnIndex,1);
y_extraFractionOfGridSpace = rem(gridRowIndex,1);

% convert the fractions of a grid space above into pixels
x_pixelsBetweenVertices = x_gridVertexToRight - x_gridVertexToLeft;
y_pixelsBetweenVertices = y_gridVertexBelow -y_gridVertexAbove;
x_pixelsToRightOfLeftVertex = x_extraFractionOfGridSpace * x_pixelsBetweenVertices;
y_pixelsBelowHigherVertex  =  y_extraFractionOfGridSpace * y_pixelsBetweenVertices;

% get pixel screen coordinates of specified input location
screen_coordinates(1) = x_gridVertexToLeft + x_pixelsToRightOfLeftVertex; % x coordinate in pixels
screen_coordinates(2) = y_gridVertexAbove + y_pixelsBelowHigherVertex; % y coordinate in pixels

