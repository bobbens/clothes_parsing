% extract contours and neighboring regions given non-max suppressed edge map
function contours = contour(nmax)

% extract contours
tic;
[skel, labels, is_v, is_e, assign, vertices, edges, ...
 v_left, v_right, e_left, e_right, c_left, c_right, ...
 edge_equiv_ids, is_compl, e_x_coords, e_y_coords] = ...
   mex_contour_sides(nmax,true);
toc;

% store pixel assignment maps
contours.skel   = skel;
contours.labels = labels;
contours.is_v   = is_v;
contours.is_e   = is_e;
contours.assign   = assign + 1;     % adjust from 0 to 1 based indexing

% store vertices and edges
contours.vertices = vertices + 1;   % adjust from 0 to 1 based indexing
contours.edges    = edges + 1;      % adjust from 0 to 1 based indexing

% store coordinates of pixels on edges (excluding endpoints)
contours.edge_x_coords = e_x_coords;
contours.edge_y_coords = e_y_coords;
for n = 1:length(e_x_coords)
   contours.edge_x_coords{n} = contours.edge_x_coords{n} + 1; % 0 -> 1 indexing 
   contours.edge_y_coords{n} = contours.edge_y_coords{n} + 1; % 0 -> 1 indexing
end

% store edge equiv ids
contours.edge_equiv_ids = edge_equiv_ids + 1;

% store completion flags
contours.is_completion = is_compl;

% get image size
im_size = size(skel);

% store vertex regions
contours.regions_v_left  = convert_cell_inds(v_left,  im_size);
contours.regions_v_right = convert_cell_inds(v_right, im_size);

% store edge regions
contours.regions_e_left  = convert_cell_inds(e_left,  im_size);
contours.regions_e_right = convert_cell_inds(e_right, im_size);

% store contour regions
contours.regions_c_left  = convert_cell_inds(c_left,  im_size);
contours.regions_c_right = convert_cell_inds(c_right, im_size);

% convert from c to matlab indices
function cell_inds = convert_cell_inds(cell_inds, im_size)
size_x = im_size(1);
size_y = im_size(2);
for n = 1:length(cell_inds)
   cell_inds{n} = convert_inds(cell_inds{n}, size_x, size_y);
end

% convert from c to matlab indices
function inds = convert_inds(inds, size_x, size_y)
x = floor(inds./size_y) + 1;
y = mod(inds,size_y) + 1;
inds = sub2ind([size_x size_y], x, y);
