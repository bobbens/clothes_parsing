function Example01()

%addpath('../x64/Debug/');
%addpath('../x64/Release/');
addpath('..');

Samples{1}.VariableCardinalities = [2 2 2 2 2 2 2 2];
Samples{1}.MachineIDs = [1 1 1 1 2 2 2 2];%ignored when using matlab!!!

Parents = {[9 10],[9 11],[10 12 13],[11 12 14],[13 15 16],[14 15 17],[16 18],[17 18]};
FactorVariables = {[1 2],[1 3],[2 4],[3 4],[3 5],[4 6],[5 6],[5 7],[6 8],[7 8]};

for k=1:8
    Samples{1}.Regions{k}.c_r = 1;
    Samples{1}.Regions{k}.VariableIndices = [k];
    Samples{1}.Regions{k}.Features{1}.r = k;
    Samples{1}.Regions{k}.Features{1}.pot = sparse([1;-1]);
    Samples{1}.Regions{k}.Loss = sparse([0;0]);
    Samples{1}.Regions{k}.Parents = Parents{k};
end
for k=1:10
    Samples{1}.Regions{k+8}.c_r = 1;
    Samples{1}.Regions{k+8}.VariableIndices = FactorVariables{k};
    Samples{1}.Regions{k+8}.Features{1}.r = 9;
    Samples{1}.Regions{k+8}.Features{1}.pot = sparse([1;-1;-1;1]);
    Samples{1}.Regions{k+8}.Loss = sparse([0;0;0;0]);
    Samples{1}.Regions{k+8}.Parents = [];
end

Samples{1}.Observation = [1 1 2 2 2 2 1 1];

Samples{2} = Samples{1};

%WriteOutput(Samples);
%PlotGraph(Samples,1);

Params.ArmijoIterations = 50;
Params.C = 1;
Params.CRFEraseMessages = 1;
Params.CRFIterations = 20;
Params.CRFMPIterations = 50;
Params.epsilon = 1.0;
Params.p = 2;
Params.ReuseMessagesForF = 0;

task = 1;
w_init = ones(9,1);

[w, Prediction] = structuredPrediction(Samples, Params, task, w_init);
%[w, Prediction] = latentStructuredPrediction(Samples, Params, task, w_init);
stop = 1;
