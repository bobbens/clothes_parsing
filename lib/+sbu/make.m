function make(varargin)
%MAKE build dependency files
%
%    sbu.ClothParser.make(...)
%
bsr.make(varargin{:});
libsvm.make(varargin{:});
liblinear.make(varargin{:});
libdai.make(varargin{:});
uci.make(varargin{:});
end
