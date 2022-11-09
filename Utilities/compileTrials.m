clear all
tLoc = 'OU_Datasets'
dLoc = dir(tLoc)
dLoc = dLoc(3:end)

structB = table();
ii = 1;
for Q = 1:length(dLoc)
    Q
    load(tLoc + "/" + dLoc(Q).name);
    for j = 1:length(mRes)
        mRes(j).TrialNumber = 1;
    end
    a = struct2table(mRes);
    structB = [structB; a];
end

clearvars -except structB tLoc
mRes = table2struct(structB);

clear structB;

save(string(tLoc) + '_compiled.mat')
done = 1