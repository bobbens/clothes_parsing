function Convert2Dot(fn_feature, fn_graph)

fid = fopen(fn_feature,'rb');
numVars = fread(fid,1,'int32');
CardMachID = fread(fid,2*numVars,'int32');
numRegions = fread(fid,1,'int32');
Regions = cell(1,numRegions);
adj = sparse(numRegions,numRegions);
node_label = cell(1,numRegions);
for r=1:numRegions
    id = fread(fid,1,'int32');
    node_label{r} = num2str(id);
    Regions{r}.c_r = fread(fid,1,'double');
    numVarsInRegion = fread(fid,1,'int32');
    Regions{r}.VariableIndices = fread(fid,numVarsInRegion,'int32')+1;
    numFeatures = fread(fid,1,'int32');
    FeatureLossInfo = fread(fid,2*numFeatures+1,'int32');
    numParents = fread(fid,1,'int32');
    Regions{r}.Parents = [];
    if(numParents>0)
        Regions{r}.Parents = fread(fid,numParents,'int32')+1;
    end
    adj(r,Regions{r}.Parents) = 1;
    adj(Regions{r}.Parents,r) = 1;
end
fclose(fid);

addpath('GraphViz2Mat1.2');
graph_to_dot(adj,'filename',fn_graph,'node_label',node_label,'directed',0)