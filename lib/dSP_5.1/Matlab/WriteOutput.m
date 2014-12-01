function WriteOutput( Path, Samples )

for x=1:numel(Samples)
    fid = fopen([Path,'Sample_',num2str(x),'.feature'],'wb');
    fprintf(fid,'MARKOV');
    fwrite(fid,10,'char');
    fprintf(fid,'%d', numel(Samples{x}.VariableCardinalities));
    fwrite(fid,10,'char');
    fprintf(fid,'%d ', Samples{x}.VariableCardinalities);
    fwrite(fid,10,'char');
    fprintf(fid,'%d ', Samples{x}.MachineIDs-1);
    for r=1:numel(Samples{x}.Regions)
        fwrite(fid,10,'char');
        fprintf(fid, '%d ', r-1);
        fprintf(fid, '%f ', Samples{x}.Regions{r}.c_r);
        fprintf(fid, '%d ', numel(Samples{x}.Regions{r}.VariableIndices));
        fprintf(fid, '%d ', Samples{x}.Regions{r}.VariableIndices-1);
        %fprintf(fid, '%d ', numel(Samples{x}.Regions{r}.Features));
        %for k=1:numel(Samples{x}.Regions{r}.Features)
        %    fprintf(fid, '%d ', Samples{x}.Regions{r}.Features{k}.r-1);
        %    fprintf(fid, '%f ', Samples{x}.Regions{r}.Features{k}.pot);
        %end
        nfeatures = size(Samples{x}.Regions{r}.pot,2);
        fprintf(fid, '%d ', nfeatures);
        for k=1:nfeatures;
            fprintf(fid, '%d ', Samples{x}.Regions{r}.r(k)-1);
            fprintf(fid, '%f ', Samples{x}.Regions{r}.pot(:,k));
        end
        fprintf(fid, '%f ', Samples{x}.Regions{r}.Loss);
        fprintf(fid, '%d ', numel(Samples{x}.Regions{r}.Parents));
        fprintf(fid, '%d ', Samples{x}.Regions{r}.Parents-1);
    end
    fclose(fid);
    
    fid = fopen([Path,'Sample_',num2str(x),'.observation'],'wb');
    fprintf(fid,'%d', numel(Samples{x}.Observation));
    fwrite(fid,10,'char');
    fprintf(fid,'%d ', Samples{x}.Observation-1);
    fclose(fid);
end
