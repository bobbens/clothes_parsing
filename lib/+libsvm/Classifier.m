classdef Classifier < handle
    %LIBSVM.CLASSIFIER abstract class for libsvm classifier
    properties(Abstract)
        model
    end
    
    properties (Constant,Hidden)
        PATH = ['tmp',filesep,'classifiers',filesep]
    end
    
    properties
        training_data_path
        %pointer to the cached training data
    end
    
    properties(Dependent,Hidden)
        training_data
        %accessor to the cached training data
    end
    
    methods(Abstract)
        train(this, y, x, varargin)
        R = predict(obj, y, x, varargin)
    end
    
    methods (Abstract,Static)
        this = loadobj(obj)
    end
    
    methods
        function this = Classifier(varargin)
            %CLASSIFIER abstract constructor
            %
            %  obj = libsvm.Classifier(y,x,'option',value,...)
            %  obj = libsvm.Classifier('option',value,...)
            
            % Get input format
            if nargin > 1 && isnumeric(varargin{2}) && ...
                    (isnumeric(varargin{1})||islogical(varargin{1}))
                % Train if data is given
                this.train(varargin{:});
            else
                % Only set options
                this.update_options(varargin{:});
            end
        end
        
        function update_options(this, varargin)
            %UPDATE_OPTIONS
            props = properties(this)';
            if isscalar(varargin)
                arg = varargin{1};
                if isstruct(arg)
                    for prop = fieldnames(arg)'
                        if any(strcmp(prop{:},props))
                            this.(prop{:}) = arg.(prop{:});
                        end
                    end
                end
            else
                for i = 1:2:length(varargin)
                    ind = strcmp(varargin{i},props);
                    if any(ind)
                        this.(props{find(ind,1)}) = varargin{i+1};
                    end
                end
            end
        end
        
        function s = struct(this)
            %struct
            s = struct;
            m = metaclass(this);
            if isprop(m,'PropertyList'), m = m.PropertyList;
            else m = [m.Properties{:}]; end
            for i = 1:numel(this)
                s(i).class = class(this);
                for prop = {m([m.Constant]==false&[m.Dependent]==false).Name}
                    s(i).(prop{:}) = this(i).(prop{:});
                end
            end
        end
        
        function obj = saveobj(this)
            %SAVEOBJ
            obj = this.to_struct;
        end
        
        function obj = copy(this)
            %COPY
            obj = feval(class(this),this.to_struct);
        end
        
        function set.training_data(this, X)
            %SET.TRAINING_DATA save or clear training data
            if isempty(X)
                if ~isempty(this.training_data_path) &&...
                        exist(this.training_data_path,'file')
                    delete(this.training_data_path);
                end
                this.training_data_path = [];
            else
                if ~exist(this.PATH,'dir'), mkdir(this.PATH); end
                if isempty(this.training_data_path)
                    this.training_data_path =...
                        [tempname(sbu.ClothParser.PATH),'.mat'];
                end
                save(this.training_data_path,'X');
            end
        end
        
        function X = get.training_data(this)
            %GET.TRAINING_DATA retrieve training data from cache
            if ~isempty(this.training_data_path) &&...
                    exist(this.training_data_path,'file')
                load(this.training_data_path,'X');
            else
                X = [];
            end
        end
    end
    
    methods(Static)
        function [y, x] = read(file_path)
            %READ read labels and vectors from libsvm format file
            [y, x] = libsvm.libsvmread(file_path);
        end
        
        function [y, x] = write(file_path, y, x)
            %WRITE write labels and vectors into libsvm format file
            if ~issparse(x), x = sparse(x); end
            libsvm.libsvmwrite(file_path, y, x);
        end
    end
    
end

