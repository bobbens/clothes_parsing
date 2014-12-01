classdef Segmentation < handle
   %SEGMENTATION segmentation structure
   properties
      map_path = '' % path to the saved map
      map         % map of superpixel indices (cached property)
      length      % size of superpixels
      labels      % labels for each superpixel
      features     % feature vectors for each superpixel
      pairs       % pairs of neighboring superpixels
      area        % area of each superpixel
   end
   
   properties (Constant)
      THRESH = 0.05 % threshold value used in Berkeley segmentation
   end
   
   methods
      %% OBJECT INITIALIZATION AND ACCESSORS
      function this = Segmentation( varargin )
         %SEGMENTATION
         %
         %   this = sbu.Segmentation(im, ...)
         %   this = sbu.Segmentation(...)
         %
         if nargin>0;
            if isnumeric(varargin{1});
               im       = varargin{1};
               varargin = varargin(2:end);
               s        = bsr.Segmentation(im);
               this.map = s.labels(this.THRESH);
            elseif isstruct(varargin{1});
               %inseg    = varargin{1};
               %this.map = inseg.map;
               %this.labels = inseg.labels;
               %varargin = varargin(2:end);
               S = varargin{1};
               props = intersect(properties(this),fieldnames(S));
               for i = 1:numel(props)
                  this.(props{i}) = S.(props{i});
               end
            end
            this.initialize(varargin{:});
         end
      end
      
      function initialize( this, varargin )
         %INITIALIZE
         if isscalar(varargin) && isstruct(varargin{1})
            S = varargin{1};
         else
            S = struct(varargin{:});
         end
         props = intersect(properties(this),fieldnames(S));
         for i = 1:numel(props)
            this.(props{i}) = S.(props{i});
         end
      end
      
      function [ s ] = struct(this)
         %STRUCT convert to a struct array
         s = struct( 'map_path',{this.map_path},...
                  'map',    cell(1,numel(this)),...
                  'labels',  {this.labels},...
                  'features',{this.features});
         % load map when a storage file does not exist
         ind = arrayfun(@(x)~exist(x.map_path,'file'),s);
         [s(ind).map] = deal(this(ind).map);
      end
      
      function obj = copy(this)
         %COPY clone an instance
         obj = sbu.Segmentation(this.struct);
      end
      
      function x = get.map(this)
         %GET.MAP lazy loading of map data
         if isempty(this.map)
            this.map = this.load_map;
         end
         x = this.map;
      end
      
      function x = get.labels(this)
         %GET.LABELS
         if isempty(this.labels)
            this.labels = zeros(this.length,1,'uint32');
         end
         x = this.labels;
      end
      
      function x = get.length(this)
         %GET.LENGTH lazy calculation of superpixel size
         if isempty(this.length) || this.length==0
            this.length = numel(unique(this.map(:)));
         end
         x = this.length;
      end
      
      function x = get.pairs(this)
         %GET.PAIRS pair indices
         if ~isempty(this.map) && isempty(this.pairs)
            this.pairs = this.pairs_from_map();
         end
         x = this.pairs;
      end
      
      function [ a ] = get.area(this)
         %GET.AREA calculate areas of each segment
         if isempty(this.area) && ~isempty(this.map)
            this.area = full(sparse(double(this.map(:)),1,1));
         end
         a = this.area;
      end
      
      %% API
      function x = mask_of(this, id)
         %MASK_OF return logical mask of the segment
         x = (this.map == id);
      end
      
      function features_from(this, fmap, aggregator)
         %FEATURES_OF segmentation
         %  fmap:      M-by-N-by-D array of features.
         %  aggregator: Aggregation function applied to set of
         %           pixelwise features. The function should take
         %           row vectors of pixel representation and return
         %           a single row vector of the segment feature.
         %           (default: @(x) [mean(x),std(x)])
         if nargin < 3, aggregator = @(x) [mean(x),std(x)]; end
         
         % format fmap into row vectors
         fmap = reshape(fmap,[size(fmap,1)*size(fmap,2),size(fmap,3)]);
         % get segment features by aggregating a set of pixel features
         ind = accumarray(this.map(:),1:numel(this.map),...
            [this.length,1],@(x){x});
         for i = 1:numel(ind)
            %m = this.mask_of(uint32(i));
            m = ind{i};
            if iscell(aggregator)
               x = fmap(m(:),:);
               f = cell2mat(cellfun(@(j) aggregator{j}(x(:,j)),...
                     1:size(fmap,3), 'UniformOutput',false));
            else
               f = aggregator(fmap(m(:),:));
            end
            if i==1 % allocate memory
               this.features = zeros(this.length,numel(f));
            end
            this.features(i,:) = f;
         end
      end
      
      function [ L ] = label_map(this)
         %LABEL_MAP return map of labels
         L = uint16(this.labels(this.map));
      end
      
      function [ ind ] = is_overlapped(this, bbox, varargin)
         %IS_OVERLAPPED return trimap given bbox
         
         thresh = 0.75;
         for i = 1:2:numel(varargin)
            switch varargin{i}
               case 'thresh', thresh = varargin{i+1};
            end
         end
         
         % get roi
         x0 = max(1,bbox(1));
         y0 = max(1,bbox(2));
         x1 = min(size(this.map,2),bbox(1)+bbox(3)-1);
         y1 = min(size(this.map,1),bbox(2)+bbox(4)-1);
         roi_map = this.map(y0:y1,x0:x1);
         
         % threshold on overlapping ratio for each
         segment_ids = unique(roi_map(:));
         ind = false(this.length,1);
         for i = 1:numel(segment_ids)
            segment_area = nnz(this.map(:)==segment_ids(i));
            roi_area = nnz(roi_map(:)==segment_ids(i));
            ind(segment_ids(i)) = roi_area / segment_area > thresh;
         end
      end
      
      %% Utility
      function [ I ] = show_superpixels(this, im, varargin)
         %SHOW_SUPERPIXELS
         if nargin == 2
            % When given an image
            edge_color = [255,0,0];
            m = this.map;
            e = (m ~= [m(2:end,:);m(end,:)]) |...
               (m ~= [m(1,:);m(1:end-1,:)]) |...
               (m ~= [m(:,2:end),m(:,end)]) |...
               (m ~= [m(:,1),m(:,1:end-1)]);
            for i = 1:size(im,3)
               I = im(:,:,i);
               I(e) = edge_color(i);
               im(:,:,i) = I;
            end
            I = im;
            h = imshow(I);
         else
            % When nothing is given
            I = label2rgb(this.map);
            h = imshow(I);
         end
      end
      
      function [ I ] = show(this, varargin)
         %SHOW show and return 2D array of labels
         
         % Options
         I = [];
         clothings = [];
         show_colorbar = true;
         ALPHA = .67;
         for i = 1:2:numel(varargin)
            switch varargin{i}
               case 'image', I = varargin{i+1};
               case 'labels', clothings = varargin{i+1};
               case 'alpha',  ALPHA = varargin{i+1};
               case 'colorbar', show_colorbar = varargin{i+1};
            end
         end
         if ischar(I), I = imread(I); end
         L = this.labels;
         
         % Create index
         labels_ = sort(unique(L(:)));
         L = arrayfun(@(l) find(labels_==l),L);
         L = L(this.map);
         %hsvmap  = hsv(numel(clothings)-1);
         %cmap = [ones(1,3);hsvmap(hsvorder,:)];
         %cmap = [ones(1,3);distinguishable_colors(numel(clothings)-1, [ALPHA ALPHA ALPHA])];
         cmap = [...
    1.0000    1.0000    1.0000; ... 1
         0         0         0;
    1.0000    0.6552    0.9655;
    1.0000         0         0;
    1.0000         0    0.7241; ... 5
    1.0000    0.8276         0;
         0         0    0.3103;
    0.3448         0         0;
         0    0.4483         0;
         0    1.0000    0.8276; ... 10
         0    0.5862    1.0000;
         0         0    1.0000;
    1.0000    0.5862    0.4138;
    0.6207    0.3448    0.9310;
    0.7586    1.0000    0.5172; ... 15
         0    1.0000         0;
    0.8966    0.0345    0.3448;
    0.1379    0.8276    0.9655;
    0.5517    0.1379    0.4138;
    0.5517    0.4828         0;
         0    0.5517    0.4483;
    0.3448    0.3448    0.5517;
    1.0000    0.8966    0.6552;
    0.9310         0    1.0000;
    0.2414    0.2069         0;
         0         0    0.6207;
    0.0345    1.0000    0.5172;
    0.8276    1.0000         0;
    0.6897    0.1724    0.0345 ];
         cmap = [...
    1.0000    1.0000    1.0000; ... 1
         0         0    1.0000;
         0    1.0000         0;
    1.0000         0         0;
    1.0000         0    0.7241;
    1.0000    0.8276         0;
         0         0    0.3103;
    0.3448         0         0;
         0    0.4483         0;
         0    0.5862    1.0000;
         0    1.0000    0.8276;
         0    0.1379    0.1379;
    1.0000    0.5862    0.4138;
    0.6207    0.3448    0.9310;
    0.7586    1.0000    0.5172;
    0.8966    0.0345    0.3448;
    0.5517    0.1379    0.4138;
    0.1379    0.8276    0.9655;
    0.5517    0.4828         0;
         0         0         0;
         0    0.5517    0.4483;
    0.3448    0.3448    0.5517;
    1.0000    0.8966    0.6552;
    0.9310         0    1.0000;
    0.2414    0.2069         0;
         0         0    0.6207;
    0.0345    1.0000    0.5172;
    0.8276    1.0000         0;
    0.6897    0.1724    0.0345;
    0.8276    0.4828    0.5517;
    0.7586    1.0000    0.8966;
    1.0000    0.4828         0;
    0.1724         0    0.1379;
         0    0.2759    0.7241;
    0.6552    0.7586    1.0000;
    0.5517    0.6897    0.4138;
    0.3103    0.7586    0.0345;
    0.1034    0.4828    0.5862;
    0.5862    0.3448    0.1724;
    0.4483    0.3448    0.3448;
    0.4138         0    0.5172;
    1.0000    0.6552    0.9655;
         0    0.7241    0.4138;
    0.6207    0.5172    0.8966;
    0.6552    0.7241         0;
    1.0000    0.4483    1.0000;
    1.0000    0.3793    0.6552;
         0    0.3103    1.0000;
         0    0.2759    0.1034;
    0.5172         0    0.8966;
    1.0000    0.8621    1.0000;
    0.9310    0.6897    0.2759;
    0.6897    0.5862    0.4138;
    1.0000    0.9655    0.4483;
    0.6897    0.5517    0.7241;
    1.0000    0.7586    0.7241;
    0.6207    0.1724    0.2759;
    0.7241         0    0.6897;
    0.3793    0.4483    0.3103];
         cmap = [...
    1.0000    1.0000    1.0000; ... 1
         0         0         0;
    1.0000    0.6552    0.9655;
    1.0000         0         0;
    1.0000         0    0.7241; ... 5
    1.0000    0.8276         0;
         0         0    0.3103;
    0.3448         0         0;
         0    0.4483         0;
         0    1.0000    0.8276; ... 10
         0    0.5862    1.0000;
         0         0    1.0000;
    1.0000    0.5862    0.4138;
    0.6207    0.3448    0.9310;
    0.7586    1.0000    0.5172; ... 15
         0    1.0000         0;
    0.8966    0.0345    0.3448;
    0.1379    0.8276    0.9655;
    0.5517    0.1379    0.4138;
    0.5517    0.4828         0;
         0    0.5517    0.4483;
    0.3448    0.3448    0.5517;
    1.0000    0.8966    0.6552;
    0.9310         0    1.0000;
    0.2414    0.2069         0;
         0         0    0.6207;
    0.0345    1.0000    0.5172;
    0.8276    1.0000         0;
    0.6897    0.1724    0.0345 ];

         %cmap = [ones(1,3);distinguishable_colors(numel(clothings)-1, [ALPHA ALPHA ALPHA])];
         %cmap
         %aoeu

         
         labelid = arrayfun( @(l) find(arrayfun( @(j)ismember(l,clothings(j).id), 1:numel(clothings) )>0), labels_ );
         cmap = cmap( labelid, : );
         
         % Render
         if isempty(I)
            imshow(L,cmap);
         else
            imshow(I);
            hold on; h = imshow(L,cmap); hold off;
            set(h,'AlphaData',ones(size(L))*ALPHA);
         end
         % add colorbar
         if ~isempty(clothings)
            names = arrayfun( ...
                  @(i) clothings( ...
                        arrayfun(@(j) ismember(labels_(i),clothings(j).id), 1:numel(clothings))), 1:numel(labels_), 'UniformOutput', false );
            ind         = cellfun(@isempty,names);
            names(ind)  = repmat({'null'},1,nnz(ind));
            names(~ind) = cellfun(@(x) x.name,names(~ind),...
                  'UniformOutput',false);
            if show_colorbar
               colorbar('TickLength',[0,0],...
                      'YTick',0:numel(labels_)-1,...
                      'YTickLabel',names,...
                      'YTickMode','manual',...
                      'YTickLabelMode','manual');
            end
         end
         if nargout > 0
            %drawnow;
            g = getframe(gca);
            I = g.cdata(1:end-1,1:end-1,:);
         end
      end
      
      function c = channel_map(this, x)
         %CHANNEL_MAP convert a map into uint8 format so that the map
         %can be saved as an image
         if nargin < 2, x = this.map; end
         msk = uint32(255);
         c = uint8(cat(3,...
            bitand(x,msk),...
            bitand(bitshift(x,-8),msk),...
            bitand(bitshift(x,-16),msk)));
      end
      
      function x = dechannel_map(this, c)
         %DECHANNEL_MAP convert an uint8 format to uint32 format
         if nargin < 2, c = imread(this.map_path); end
         c = uint32(c);
         % Dechannel on load
         x = bitor(bitor(c(:,:,1),...
            bitshift(c(:,:,2),8)),...
            bitshift(c(:,:,3),16));
      end
      
      function save_map(this, file_path)
         %SAVE_MAP save a map file to a specified path
         if nargin < 2, file_path = this.map_path; end
         imwrite(this.channel_map,file_path);
         this.map_path = file_path;
      end
      
      function m = load_map(this, file_path)
         %LOAD_MAP load a map file from a specified path
         if nargin > 1, this.map_path = file_path; end
         m = [];
         if ischar(this.map_path) && exist(this.map_path,'file')
            m = this.dechannel_map;
         end
      end


      function [indices, boxes, neighbours] = blobs( this )
         indices = double(this.map);
         boxes = zeros( this.length, 4 );
         for x=1:size(this.map,1);
            for y=1:size(this.map,2);
               i = this.map(x,y);
               if boxes(i,1) > x; boxes(i,1) = x; end
               if boxes(i,2) > y; boxes(i,2) = y; end
               if boxes(i,3) < x; boxes(i,3) = x; end
               if boxes(i,4) < y; boxes(i,4) = y; end
            end
         end
         neighbours = diag(ones(this.length,1));
         for i=1:size(this.pairs,1);
            m = this.pairs(i,1);
            n = this.pairs(i,2);
            neighbours(m,n) = 1;
            neighbours(n,m) = 1;
         end
      end
   end
   
   %% Internal
   methods (Access = private)
      function [ s ] = pairs_from_map( this )
         %PAIRS_FROM_MAP find neighboring pairs of segment ids
         
         % get boundary index
         m = this.map;
         shift.u = (m ~= [m(2:end,:);m(end,:)]);
         shift.d = (m ~= [m(1,:);m(1:end-1,:)]);
         shift.l = (m ~= [m(:,2:end),m(:,end)]);
         shift.r = (m ~= [m(:,1),m(:,1:end-1)]);
         % get neighboring segment ids in row vectors
         s = [m(shift.u),m(shift.d);m(shift.l),m(shift.r)];
         s = unique(s,'rows');
         s = unique([s;fliplr(s)],'rows');
         s(s(:,1)>s(:,2),:) = [];
      end
   end
   
   %% Helpers
   methods (Hidden)
      function obj = saveobj(this)
         %SAVEOBJ callback on save
         obj = this.struct;
      end
   end
   
   methods (Static)
      function [ P ] = color_palette(L, clothings)
         %COLOR_PALETTE return color palette
         labels_ = unique(L(:));
         cmap = [ones(1,3);hsv(numel(labels_)-1)];
         P = struct('r',num2cell(cmap(:,1)),...
                  'g',num2cell(cmap(:,2)),...
                  'b',num2cell(cmap(:,3)));
         if nargin > 1
            names = arrayfun(@(i)clothings(ismember(labels_i(i), [clothings.id])),...
               1:numel(labels_),'UniformOutput',false);
            ind = cellfun(@isempty,names);
            names(ind)  = repmat({'null'},1,nnz(ind));
            names(~ind) = cellfun(@(x) x.name,names(~ind),...
               'UniformOutput',false);
            if isempty(names), names = {'null'}; end
            [P.name] = deal(names{:});
         end
      end
   end
   
   methods (Hidden, Static)
      function [ this ] = loadobj(obj)
         %LOADOBJ callback on load
         this = sbu.Segmentation(obj);
      end
   end
   
end

