classdef Quantizer < handle
    %QUANTIZER abstract quantizer class
    
    properties(Abstract)
        model
        % This is where the learned model is stored
        % This property can be any matlab data
    end
    
    methods
        function [ this ] = Quantizer( varargin )
            %QUANTIZER constructor used in subclasses
            props = properties(this);
            if isscalar(varargin)
                s = varargin{1};
                if isstruct(s)
                    fields = fieldnames(s);
                    for i = 1:numel(fields)
                        if any(strcmp(fields{i},props))
                            this.(fields{i}) = s.(fields{i});
                        end
                    end
                end
            else
                for i = 1:2:length(varargin)
                    ind = strcmp(varargin{i},props);
                    if any(ind), this.(props{ind}) = varargin{i+1}; end
                end
            end
        end
        
        function [ s ] = to_struct( this )
            %TO_STRUCT
            s = struct;
            props = properties(this);
            for i = 1:numel(props)
                s.(props{i}) = this.(props{i});
            end
        end
        
        function [ obj ] = copy(this)
            %COPY
            obj = improc.Quantizer(this.to_struct);
        end
    end
    
    methods(Abstract)
        [ idx ] = quantize( this, vectors, varargin )
        % Give indices to data samples using a learned model
        % The function should take row vectors and return a vector of
        % quantized results
        learn( this, vectors, varargin )
        % Learn a model from data samples
        % This function should create a new 'model' property
    end
    
end

