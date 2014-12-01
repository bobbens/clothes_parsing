function startup
   %STARTUP set up path information on matlab startup
   root = fileparts(mfilename('fullpath'));
   addpath(fullfile(root, 'lib'));
   for d = dir(fullfile(root, 'lib'))'
      if d.isdir && d.name(1) ~= '.' && d.name(1) ~= '+'
         addpath(fullfile(root, 'lib', d.name));
      end
   end
   addpath(fullfile(root,'tasks'));
   addpath(fullfile(root,'lib/SelectiveSearchCodeIJCV/Dependencies'));
end
