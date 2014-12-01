% interactively display contours and neighboring regions
function disp_contours(contours, im)

% get image size
im_size = size(contours.skel);

% get vertex and edge indices
v_inds = find(contours.is_v);
e_inds = find(contours.is_e);

% initialize region indices
r_inds_left  = [];
r_inds_right = [];
r_inds_left_ext  = [];
r_inds_right_ext = [];

% select and show edge regions
while (1)
   % initialize pixel color maps
   r = zeros(im_size);
   g = zeros(im_size);
   b = zeros(im_size);
   % highlight extended regions
   b(r_inds_left_ext)  = 0.50;
   r(r_inds_right_ext) = 0.50;
   % highlight regions
   b(r_inds_left)  = 1;
   r(r_inds_right) = 1;
   % highlight edges
   g(e_inds) = 1;
   % highlight vertices
   r(v_inds) = 0.75;
   g(v_inds) = 0.25;
   % assemble color map
   pixel_type(:,:,1) = r;
   pixel_type(:,:,2) = g;
   pixel_type(:,:,3) = b;
   % display
   if (nargin > 1)
      % display vertices, edges, and regions
      figure(1); imagesc(pixel_type.*0.5 + im.*0.5); axis image;
      % display image clipped to regions
      region_inds = unique( ...
         [r_inds_left; r_inds_right; r_inds_left_ext; r_inds_right_ext] ...
      );
      if (~isempty(region_inds))
         [xs ys] = ind2sub(im_size, region_inds);
         inds_r = sub2ind([im_size 3],xs,ys,ones(size(xs)));
         inds_g = sub2ind([im_size 3],xs,ys,2*ones(size(xs)));
         inds_b = sub2ind([im_size 3],xs,ys,3*ones(size(xs)));
         im_clip = zeros(size(im));
         im_clip(inds_r) = im(inds_r);
         im_clip(inds_g) = im(inds_g);
         im_clip(inds_b) = im(inds_b);
         x_min = min(xs); x_max = max(xs);
         y_min = min(ys); y_max = max(ys);
         im_clip = im_clip(x_min:x_max,y_min:y_max,:);
         figure(2); imagesc(im_clip); axis image;
      else
         figure(2); imagesc(zeros(size(im))); axis image;
      end
   else
      % display vertices, edges, and regions
      figure(1); imagesc(pixel_type); axis image;
   end
   % get click
   figure(1);
   [y, x, button] = ginput(1);
   if (button == 1)
      % select vertex/edge clicked by user
      x = round(x);
      y = round(y);
      if ((x > 0) && (x < im_size(1)) && (y > 0) && (y < im_size(2)))
         if (contours.is_v(x,y) == 1)
            % get vertex id
            v_id = contours.assign(x,y);
            % reset extended regions
            r_inds_left_ext  = [];
            r_inds_right_ext = [];
            % set left/right regions of vertex
            r_inds_left  = contours.regions_v_left{v_id};
            r_inds_right = contours.regions_v_right{v_id};
         elseif (contours.is_e(x,y) == 1)
            % get edge and contour ids
            e_id     = contours.assign(x,y);
            equiv_id = contours.edge_equiv_ids(e_id);
            % set extended contour regions
            r_inds_left_ext  = contours.regions_c_left{equiv_id};
            r_inds_right_ext = contours.regions_c_right{equiv_id};
            % set edge regions
            r_inds_left  = contours.regions_e_left{e_id};
            r_inds_right = contours.regions_e_right{e_id};
         end
      end
   elseif (button == 3)
      % reset extended regions
      r_inds_left_ext  = [];
      r_inds_right_ext = [];
      % reset regions 
      r_inds_left  = []; 
      r_inds_right = [];
   else
      break;
   end
end
