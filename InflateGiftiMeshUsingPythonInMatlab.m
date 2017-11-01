function imesh = InflateGiftiMeshUsingPythonInMatlab(gifti_file)
% Inflate a gifti mesh from within MATLAB using python
%
% - Requirements:
%   - matlab 'gifti' constructor
%   - python or anaconda installed
%   - my PyBP mesh tools
% AS


if nargin < 1
    gifti_file = 'spm.surf.gii';
end

mesh = gifti(gifti_file);


% PyBP [& conda] in MATLAB
% pyversion /Users/Alex/anaconda/bin/python

% PyBP location
modpath = '/Users/Alex/code/PyBP';
P       = py.sys.path;
if count(P,modpath) == 0
    insert(P,int32(0),modpath);
end

% Import the pybp tools
pybp = py.importlib.import_module('PyBP');
vf   = pybp.GetMesh(gifti_file);
infl = pybp.inflate(vf{1},vf{2});

% Rearrange returned data to matlab double
data = double(py.array.array('d',py.numpy.nditer(infl))); 
shp  = infl.shape;
shp  = shp.cell; 
data = reshape(data,[shp{2} shp{1}])'; 

imesh = mesh;
imesh.vertices = data;