function Samples = ReadOutput(fn)

fid = fopen(fn,'rb');
task = fread(fid,1,'int32');

if(task<2)
    numFeat = fread(fid,1,'int32');
    w = fread(fid,numFeat,'double');
end
if(task>0)
    numSamples = fread(fid,1,'int32');
    Samples = cell(1,numSamples);
    for x=1:numSamples
        numRegions = fread(fid,1,'int32');
        Samples{x}.Regions = cell(1,numRegions);
        for r=1:numRegions
            numBel = fread(fid,1,'int32');
            Samples{x}.Regions{r}.bel = fread(fid,numBel,'double');
        end
    end
end
fclose(fid);