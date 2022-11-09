%% HEADER

% A002_Generate_OU_Dataset.m
% From: An emergent temporal basis set robustly supports cerebellar time-series learning
% Authors: Jesse I. Gilmer, Michael A. Farries, Zachary Kilpatrick, Ioannis Delis, Abigail L. Person
% 2022, BiorXiv, https://www.biorxiv.org/content/10.1101/2022.01.06.475265v1
%
% This script was written by JG, with help from MAF.
%
% INPUTS: N/A.
% OUTPUTS: This script creates a series of datasets. They can be combined
% into a single file for analysis by Compile_Datasets.m


%% BODY

% clear worksapce
clear all
close all

% setup save locale
save_dir = 'OU_Datasets/';
mkdir(save_dir)

% Cycle through thresholds:
for A = [-2:.1:2]
    clearvars -except save_dir A
    warning('off')

    % create the output structure.
    mRes = struct();
    
    % Network params:
    nMF = 50;   % number of mossy fibers. Default = 50.
    nGC = 300; % number of granule cells. Default = 3000.
    conv = 4;   % convergence ratio, mf:gc. Default = 4.
    
    % MF params:
    mu = 4;    % Mean rate. Default = 4, but this value is normalized. 
    sigma = 1; % Rate standard deviation. Default = 1, but this value is normalized. 
    T = 100;   % OU Tau, the decay rate. Default = 100.
    
    % Integrate into the options storage dictionary.
    options.mu = mu;
    options.sigma = sigma;
    
    % Set the learn rate and trials.
    options.lrnRate = 1e-5; %Learning rate. This will change below.
    options.nTrials = 100; %Number of trials to learn over. Default = 1000.
    
    % Set the length of the target function:
    nPts = 1000; %Default 1000
    
    % Set the number of independent simulations to run:
    ntr = 50; % Default = 50.
    
    % i is the external itterator.
    i = 1;

    close all
    for q = 1:ntr
        [q/ntr A]
        %% Project params:
        alpha = A;
        mRes(i).options = options;
        options.tau = T;
        mRes(i).alpha = alpha;
        mRes(i).conv = conv;
        mRes(i).nGC = nGC;
        mRes(i).nMF = nMF;
        
        
        %% Make input signal:
        input = MakeSignal(nPts,nMF,1,options);
        nPts = size(input,2);
        mRes(i).nPts = nPts;
        
        %% Make target signal:
        
        options2 = options;
        options2.tau = 10;
        target = MakeSignal(nPts,1,1,options2);
        % target = rand(size(target));
        target = (target-min(target))./(max(target)-min(target));
        
        %% Make GCs:
        W1 = zeros(nGC, nMF);
        for k = 1:nGC
            sel = randsample(nMF, conv);
            W1(k,sel) = 1/conv;
        end
        mu = mean(mean(input));
        sigma = mean(std(input,[],2));
        thr = mu + alpha * sigma;
        
        gc = W1 * input;
        gc = gc - thr;
        gc(gc < 0) = 0;
        % gc = gc > 0;
        [~,loss,~,~] = calcWords(gc,0);
        
        mRes(i).gcPopSparse = sum(sum(gc' > 0)<=0)/nGC;
        mRes(i).gcBinarySparse = mean(sum(gc' > 0)/nPts);
        mRes(i).gcBinaryVariance = sum(var(gc' > 0))/nGC;
        
        
        CC = corr(gc');
        diagC = tril(ones(size(CC)),-1);
        mRes(i).gcCorr = nanmean(CC(diagC == 1));
        
        options.lrnRate = 1e-4;
        [results] = RunFilter(input,target,gc,options);
        
        mRes(i).mse_all = results.mse;
        mRes(i).mse = results.mse(end);
        mRes(i).gcopt = options;
        
        try
            ft = fittype('a1*exp(-x/tau1) + c');
            x = [1:length(results.mse(2:end))]';
            y = results.mse(2:end)';
            y = (y-min(y))./(max(y)-min(y)+.0001);
            y(isinf(y)) = 0;
            f = fit(x, y, ft, 'StartPoint',  [15 1 100], 'Lower', [1 0.01 1]);
            mRes(i).gcLearntau = f.tau1;
        catch
            mRes(i).gcLearntau = NaN;
        end
        
        options.lrnRate = 1e-6;
        [results] = RunFilter(input,target,input,options);
        mRes(i).MF_mse_all = results.mse;
        mRes(i).MF_mse = results.mse(end);
        mRes(i).mfopt = options;
        
        try
            ft = fittype('a1*exp(-x/tau1) + c');
            x = [1:length(results.mse(2:end))]';
            y = results.mse(2:end)';
            y = (y-min(y))./(max(y)-min(y)+.0001);
            y(isinf(y)) = 0;
            f = fit(x, y, ft, 'StartPoint',  [15 1 100], 'Lower', [1 0.01 1]);
            mRes(i).mfLearntau = f.tau1;
        catch
            mRes(i).mfLearntau = NaN;
        end
        
        % Perform the intermediate analyses:
        % tau of autocorr
        [mRes(i).gcTau] = GetSignalTau(gc');
        [mRes(i).inTau] = GetSignalTau(input');
        
        % signal dimensionality
        mRes(i).gcDim = sum(eig(cov(gc')))^2/sum(eig(cov(gc')).^2);
        mRes(i).mfDim = sum(eig(cov(input')))^2/sum(eig(cov(input')).^2);
        
        % signal pop correlations
        CC = corr(gc');
        diagC = tril(ones(size(CC)),-1);
        [mRes(i).gcCorr]  = nanmean(CC(diagC == 1));
        
        CC = corr(input');
        diagC = tril(ones(size(CC)),-1);
        [mRes(i).mfCorr]  = nanmean(CC(diagC == 1));
        
        % TS and related metrics
        [mRes(i).gcTS_rel,mRes(i).gcLoss,mRes(i).gcTS,mRes(i).gcTS_weak] = calcWords(gc,0);
        [mRes(i).mfTS_rel,mRes(i).mfLoss,mRes(i).mfTS,mRes(i).mfTS_weak] = calcWords(input,0);
        
        % Variance
        mRes(i).gcVar = sum(nanvar(gc'));
        mRes(i).mfVar = sum(nanvar(input'));
        
        % Explained variance of input
        % Put MFs and GCs in the right form:
        ssres = 0;
        sstot = 0;
        x = input';
        y = gc';
        
        B = [];
        for k=1:nMF
            b = regress(x(:,k),y);
            B = [B, b];
        end
        
        % compute the explained variance across all these realizations
        xT = x';
        sstot = sstot+nMF*sum(var(xT));
        ssres = ssres+sum(sum((x-y*B).^2));
        
        mRes(i).explainedVar = 1-ssres/sstot;
        
        [coeff,score,latent,tsquared,explained] = pca(gc);
        
        mRes(i).PCAind100 = sum(explained>=(100/(nGC)))/nGC;
        mRes(i).PCAind1 = sum(explained>=(1/(nGC)))/nGC;
        mRes(i).PCAind = sum(explained>=0)/size(explained,1);
        mRes(i).PCAvar = var(explained);
        mRes(i).PCAmean = mean(explained);
        
        i = i + 1;
    end
    
    % Save the output:
    save(save_dir + "OU_ct_" + (ntr) + "_threshold_" + A + "s01_" + randi([1 1000],1,1) + "s02_" + randi([1 1000],1,1) + ".mat")
    
    % Restart loop over threshold:
end


%% Next step...
% Call 