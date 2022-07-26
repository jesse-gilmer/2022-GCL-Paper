clear all
close all
for nGC = [10 50 100 300 600 1000:500:5000]
    for A = 0
        clearvars -except nGC A
        warning('off')
        T = 100;
        
        save_dir = 'nGC_Test/';
        mkdir(save_dir)
        
        mRes = struct();
        
        nMF = 50; %Default 50
%         nGC = 3000; %Default 3000
        conv = 4;
        
        mu = 4;
        sigma = 1;
        
        options.mu = mu;
        options.sigma = sigma;
        options.lrnRate = 1e-5;
        options.nTrials = 500; %Default 1000
        
        nPts = 1000; %Default 1000
        
        ntr = 30;
        i = 1;
        
        
        %     for A = -4:.25:2.5
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
            
            % figure(1); hold on
            % plot(results.mse)
            %
            % figure(2); hold off
            % plot(results.W*gc)
            % hold on
            % plot(target,'k')
            % title(alpha)
            
            %% Do intermediate analyses:
            % tau of autocorr
            [mRes(i).gcTau] = GetSignalTau(gc');
            [mRes(i).inTau] = GetSignalTau(input');
            
            % figure(1)
            % subplot(1,2,1); hold on
            % scatter(A,mRes(i).gcTau,'r')
            % scatter(A,mRes(i).inTau,'k')
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
    
    %     [mRes(i-1).inTau mRes(i-1).gcTau]
    %     [T A]
    save(save_dir + "nGC_" + nGC + "Indep_PARAMS_t_" + (ntr) + "trials_threshold_" + A + "seed_" + randi([1 1000000],1,1) + "seed2_" + randi([1 1000000],1,1) + ".mat")
    SAVED = "JUST SAVED, YAY."
    end
end
