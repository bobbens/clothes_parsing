function WriteOutputBinary(Samples)

for x=1:numel(Samples)
    fid = fopen(['Sample_',num2str(x),'.feature'],'wb');
    fwrite(fid,numel(Samples{x}.VariableCardinalities),'int32');
    fwrite(fid,Samples{x}.VariableCardinalities,'int32');
    fwrite(fid,Samples{x}.MachineIDs-1,'int32');
    fwrite(fid,numel(Samples{x}.Regions),'int32');
    for r=1:numel(Samples{x}.Regions)
        fwrite(fid,r-1,'int32');
        fwrite(fid,Samples{x}.Regions{r}.c_r,'double');
        
        info = zeros(1,1+numel(Samples{x}.Regions{r}.VariableIndices)+1+2*numel(Samples{x}.Regions{r}.Features)+1+1+numel(Samples{x}.Regions{r}.Parents));
        pos = 1;
        info(pos) = numel(Samples{x}.Regions{r}.VariableIndices);
        pos = pos + 1;
        info(pos:pos+numel(Samples{x}.Regions{r}.VariableIndices)-1) = Samples{x}.Regions{r}.VariableIndices-1;
        pos = pos + numel(Samples{x}.Regions{r}.VariableIndices);
        info(pos) = numel(Samples{x}.Regions{r}.Features);
        pos = pos + 1;
        for k=1:numel(Samples{x}.Regions{r}.Features)
            info(pos) = Samples{x}.Regions{r}.Features{k}.r-1;
            info(pos+1) = Samples{x}.Regions{r}.Features{k}.potIX;
            pos = pos + 2;
        end
        info(pos) = Samples{x}.Regions{r}.LossIX;
        pos = pos + 1;
        info(pos) = numel(Samples{x}.Regions{r}.Parents);
        pos = pos + 1;
        info(pos:pos+numel(Samples{x}.Regions{r}.Parents)-1) = Samples{x}.Regions{r}.Parents-1;
        
        fwrite(fid,info,'int32');
    end
    fwrite(fid,numel(Samples{x}.Potentials),'int32');
    for p=1:numel(Samples{x}.Potentials)
        fwrite(fid,Samples{x}.Potentials{p}.IX,'int32');
        fwrite(fid,numel(Samples{x}.Potentials{p}.pot),'int32');
        fwrite(fid, Samples{x}.Potentials{p}.pot, 'double');
    end
    fclose(fid);
    
    fid = fopen(['Sample_',num2str(x),'.observation'],'wb');
    fprintf(fid,'%d', numel(Samples{x}.Observation));
    fwrite(fid,10,'char');
    fprintf(fid,'%d ', Samples{x}.Observation-1);
    fclose(fid);
end