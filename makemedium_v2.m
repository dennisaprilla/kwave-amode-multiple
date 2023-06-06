function medium_logic = makemedium_v2(theta, medium_gridsize, dx, x_pos)

% % Create rotation matrix
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
p = R(:,1);
n = R(:,2);

Nx = medium_gridsize(1);
Ny = medium_gridsize(2);

row = -Nx/2:Nx/2-1;
col = -Ny/2:Ny/2-1;
[Row, Col] = meshgrid(row, col);
% Coordinate system representing matrix and graph are different. In matrix
% row is x, while in graph row is y. So to avoid confusion, i use explicit
% naming, XY_matrix (that i don't use) and XY_graph
XY_matrix = [Row(:), Col(:)];
XY_graph  = [Col(:), Row(:)];

% While line can be represented with equation y=mx+c, here i use
% point-normal representation. Let n=(nx,ny) be a normal vector, 
% originated from a point o=(x0, y0) that lie in the line, a set of points 
% X=(x,y) that represents the line satisfies the following equation: 
% <n, (x-x0, y-y0)>, with <,> denoted as dot product. If the dot product
% results 0, those points lie in the line, if >0 above the line, if <0
% below the line

% Dot product between a point in grid (XY_graph) with normal (which will be 
% rotated from the parameter theta). 
bias_grid = round(x_pos/dx);
XY_logic = ( n'*XY_graph') > 0 + bias_grid;

% XY_logic is using graph coordinate system (row is y), we need to return
% to matrix coordinate system (row is x). First, after we reshape our 
% XY_logic to original matrix size, we transpose it. But since XY_logic
% starts from -x, we need to flip it so it starts with +x. Uncomment the
% line below to understrand what I meant.
XY_graph_logic = [XY_graph, XY_logic'];
medium_logic = flip(reshape(XY_logic, size(Row))', 1);


end

