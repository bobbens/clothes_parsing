function output = imencode(input, fmt, varargin)
%IMENCODE Encode input image using image compression algorithms.
%
%    output = imencode(input)
%    output = imencode(input, fmt)
%
% IMENCODE compresses an image input INPUT using specified format FMT.
% INPUT is M-by-N-by-d numeric array of image data. FMT is a name of image
% file extension that is recognized by IMFORMATS function, such as 'jpg' or
% 'png'. When FMT is omitted, 'png' is used as default.
%
% See also imdecode imformats imwrite

if nargin < 2, fmt = 'png'; end

tempfile = sprintf('%s.%s', tempname, fmt);
try
    imwrite(input, tempfile, fmt, varargin{:});
    fid = fopen(tempfile, 'r');
    output = fread(fid, inf, 'uint8=>uint8')';
    fclose(fid);
    delete(tempfile);
catch e
    if exist(tempfile, 'file')
        delete(tempfile);
    end
    rethrow(e);
end

end