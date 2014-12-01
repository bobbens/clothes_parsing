% extract contours and neighboring regions given non-max suppressed edge map
function contour_export(filename, nmax)

% extract contours
[skel, labels, is_v, is_e, assign, vertices, edges, e_x, e_y] = mex_contour(nmax);
n_vertices = size(vertices,1)
n_edges = size(edges,1)

% flip x,y
vertices = vertices(:,2:-1:1);

% open file 
f = fopen(filename,'w');

% write junctions
fprintf(f, '%d\n', n_vertices);
fprintf(f, '%d %d\n', vertices');

% write edges
fprintf(f, '%d\n', n_edges);
fprintf(f, '%d %d\n', edges');

% write line for each edge
for n = 1:n_edges
   coords = [e_y{n} e_x{n}];
   fprintf(f, '%d ', coords');
   fprintf(f, '-5 -5\n');
end

% write regions
fprintf(f, '0');
