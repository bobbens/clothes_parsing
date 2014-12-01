function PlotGraph(Samples, x, fn_graph)

if(nargin<3)
    fn_graph = '';
end

numRegions = numel(Samples{x}.Regions);
numVarsP1 = numel(unique(Samples{x}.MachineIDs))+1;
CMap = jet(numVarsP1);

if(numRegions>500 && nargin<3)
    return;
end

adj = zeros(numRegions);
nodeDesc = cell(1,numRegions);
nodeColors = zeros(numRegions,3);
for r=1:numRegions
    for p=1:numel(Samples{x}.Regions{r}.Parents)
        pix = Samples{x}.Regions{r}.Parents(p);
        adj(r,pix) = 1;
        adj(pix,r) = 1;
    end
    tmp = '';
    for k=1:numel(Samples{x}.Regions{r}.Features);
        tmp = sprintf('%s %d', tmp, Samples{x}.Regions{r}.Features{k}.r);
    end
    nodeDesc{r} = sprintf('c_r: %f\nFeatures: %s', Samples{x}.Regions{r}.c_r, tmp);
    varIX = Samples{x}.Regions{r}.VariableIndices;
    if(numel(varIX)==1)
        nodeColors(r,:) = CMap(Samples{x}.MachineIDs(varIX),:);
    else
        nodeColors(r,:) = CMap(end,:);
    end
end

if(nargin<3)
    g = graphViz4Matlab('-adjMat',adj,'-undirected',true,'-nodeDescriptions',nodeDesc,'-nodeColors',nodeColors);
else
    graph_to_dot(adj,'filename',fn_graph,'node_label',nodeDesc,'directed',0);
end