classdef Fashionista < sbu.Photo
    %Fashionista Fashionista dataset class
    %
    %    photos = sbu.Fashionista.load
    %
    % See also sbu.Photo
    %
    
    properties (Hidden, Constant)
        % default path to resource files
        PATH = fullfile('data','Fashionista');
        INDEX_PATH = fullfile(sbu.Fashionista.PATH,'index.mat');
        RESOURCES_PATH = fullfile(sbu.Fashionista.PATH,'resources');
        IMAGE_NAME = '%06d.jpg';
        MAP_NAME = '%06d_map.png';
    end
    
    methods
        function this = Fashionista(varargin)
            %FASHIONISTA constructor
            this@sbu.Photo(varargin{:});
        end
        
        function p = default_image_path(this)
            %DEFAULT_IMAGE_PATH
            p = fullfile(this.RESOURCES_PATH,sprintf(this.IMAGE_NAME,this.id));
        end
        
        function p = default_map_path(this)
            %DEFAULT_IMAGE_PATH
            p = fullfile(this.RESOURCES_PATH,sprintf(this.MAP_NAME,this.id));
        end
        
        function save_resources(this)
            %SAVE_RESOURCES save resources to the default path
            if ~exist(sbu.Fashionista.RESOURCES_PATH,'dir')
                mkdir(sbu.Fashionista.RESOURCES_PATH);
            end
            for i = 1:numel(this)
                im = this(i).image;
                this(i).image_path = this(i).default_image_path;
                imwrite(im, this(i).image_path);
                this(i).segmentation.save_map(this(i).default_map_path);
            end
        end
        
        function save(this)
            %SAVE save the whole dataset
            photos = this;
            photos.save_resources;
            save(sbu.Fashionista.INDEX_PATH,'photos');
        end
    end
    
    methods (Static)
        function photos = load_v2
            %LOAD_PHOTOS load the whole dataset, version 2
            S      = load(sbu.Fashionista.INDEX_PATH,'photos');
            photos = S.photos;
            for i=1:numel(photos);
               photos(i).segmentation.map_path = ...
                     [sbu.Fashionista.RESOURCES_PATH sprintf('/%06d_map.png',photos(i).id)];
            end
        end
        function photos = load
            %LOAD_PHOTOS load the entire dataset, version 3.1
            S        = load( 'data/fashionista.mat' );
            clothes  = sbu.Fashionista.clothings;
            samples  = S.samples;
            photos   = sbu.Photo();
            for i=1:numel(samples);
               % Load image
               im = imdecode( samples(i).image, 'jpg' );
               im = uint8( im );

               % Set up pose
               pose  = sbu.Pose( samples(i).poses );

               % Segmentation
               seg = samples(i).segments;
               l_  = zeros( size(seg.labels) );
               for j=1:numel(clothes);
                  l_( strcmp( clothes(j).name, seg.labels ) ) = clothes(j).id;
               end
               %segstr = struct( 'labels', {seg.labels} );
               segstr = struct( 'labels', l_', ...
                     'map', imdecode(seg.segmentation,'png') );
               segmentation = sbu.Segmentation( segstr );

               % Create image
               photos(i) = sbu.Photo( im, 'id', samples(i).id, ...
                     'pose', pose, ...
                     'segmentation', segmentation );
            end

            % Augment the data
            auglist = { 'index', 'user', 'title', 'date', 'description', 'votes', 'comments_count', 'bookmarks_count', 'photo_uri' };
            for i=1:numel(samples);
               for j=1:numel(auglist);
                  aug = auglist{j};
                  photos(i).metadata = setfield(photos(i).metadata, aug, getfield(samples(i), aug));
               end
            end
        end
        
        function C = clothings(C)
            %CLOTHINGS struct array of clothings
            % 
            % ## Input
            % * __C__ struct array to replace id/name table in the model. S
            %     must have `id` and `name` field.
            %
            % ## Output
            % * __C__ struct array of id/name table in the model.
            %
            persistent S;
            if isempty(S)
                S = load(fullfile(fileparts(mfilename('fullpath')),...
                    'clothings.mat'));
                S = S.clothings;
            end
            if nargin == 1 && isstruct(C)
                S = C;
            end
            C = S;
        end
        
        function names = clothing_names(ids)
            %CLOTHING_NAMES lookup clothing names/ids for each id/name
            %
            %    names = sbu.Fashionista.clothing_names
            %    names = sbu.Fashionista.clothing_names(ids)
            %    ids = sbu.Fashionista.clothing_names(names)
            %
            % ## Input
            % * __ids__ array of clothing ids to lookup
            % * __names__ cell array of clothing names to lookup
            %
            % ## Output
            % * __names__ cell array of clothing names found
            % * __ids__ array of clothing ids found
            %
            S = sbu.Fashionista.clothings;
            S_names = {S.name};
            S_ids = [S.id];
            
            if nargin<1, ids = S_ids; end
            if iscellstr(ids)
                names = cellfun(@(name)S_ids(strcmp(name,S_names)),ids);
            else
                names = arrayfun(@(id)S_names(id==S_ids),ids);
            end
        end
        
        function domain = validate_domain(domain)
            %VALIDATE_DOMAIN check the validity of the clothing names
            
            % Get all clothings if empty
            all_clothing_names = sbu.Fashionista.clothing_names;
            default_labels = {'null','skin','hair'};
            if nargin < 1 || isempty(domain)
                domain = all_clothing_names;
            end
            % Get id representation of names
            if iscellstr(domain)
                % Ensure default_labels
                ind = cellfun(@(x)~any(strcmp(x,domain)),default_labels);
                domain = [domain,default_labels(ind)];
                
                domain_ = intersect(domain,all_clothing_names);
                if numel(domain_) ~= numel(domain)
                    unknown = setdiff(domain,all_clothing_names);
                    unknown = sprintf('%s ',unknown{:});
                    warning('sbu:ClothParser:validateDomain',...
                        'Unknown clothing names removed: %s',unknown);
                end
                domain = sbu.Fashionista.clothing_names(domain_);
            end
            domain = sort(domain(:))';
        end
    end
    
end

