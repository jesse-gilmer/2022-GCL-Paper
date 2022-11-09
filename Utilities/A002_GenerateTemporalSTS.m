mkdir('Temporal_STS2')

clear all
close all

for gw = [2 3 4 5 10 25]
    for p = [1 0.001 0.005 .001 0.005 .01 .05 0.075 .1:.1:.9]
        clearvars -except gw p
        G = gw;
        mRes = struct();
        
        nMF = 50;
        nGC = 1000;
        conv = 4;
        
        mu = 4;
        sigma = 1;
        
        options.mu = mu;
        options.sigma = sigma;
        options.lrnRate = 3e-3;
        options.nTrials = 2000;
        
        nPts = nGC;
        
        ntr = 20;
        i = 1;
        for q = 1:ntr
            %% Project params:
            %% Project params:
            alpha = 0;
            mRes(i).options = options;
            options.tau = 10;
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
            target = rand(size(target));
            % target = (target-min(target))./(max(target)-min(target));
            
            %% Make GCs:
            gc = zeros(nPts,nPts);
            
            w = gausswin(G);
            w = w/sum(w);
            
            for k = 1:nPts
                gc(k,k) = 1;
                if (rand(1) <= p)
                    gc(k,:) = filtfilt(w,1,gc(k,:));
                end
            end
            
            gc(gc>0) = 1;
            
            subplot(2,2,1)
            imshow(gc(1:100,1:100))
            
            
            CC = corr(gc');
            diagC = tril(ones(size(CC)),-1);
            mRes(i).gcCorr = nanmean(CC(diagC == 1));
            
            [results] = RunFilter(input,target,gc,options);
            mRes(i).mse = results.mse(end);
            subplot(2,2,2)
            hold on
            plot(results.mse)
            
            ft = fittype('a1*exp(-x/tau1) + c');
            x = [1:length(results.mse(2:end))]';
            y = results.mse(2:end)';
            y = (y-min(y))./(max(y)-min(y));
            y(isinf(y)) = 0;
            f = fit(x, y, ft, 'StartPoint',  [15 1 100], 'Lower', [1 0.01 1]);
            mRes(i).gcLearntau = f.tau1;
            
            %% Do intermediate analyses:
            % tau of autocorr
            [mRes(i).gcTau] = GetSignalTau(gc');
            [mRes(i).inTau] = GetSignalTau(input');
            
            % figure(1)
            % subplot(1,2,1); hold on
            % scatter(A,mRes(i).gcTau,'r')
            % scatter(A,mRes(i).inTau,'k')
            % signal dimensionality
            mRes(i).gcDim = sum(eig(cov(gc)))^2/sum(eig(cov(gc)).^2);
            mRes(i).mfDim = sum(eig(cov(input)))^2/sum(eig(cov(input)).^2);
            
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
            subplot(2,2,3)
            hold on
            scatter(mRes(i).gcTS_rel,mRes(i).mse,5,'k');
            xlabel('STS')
            ylabel('MSE')
            subplot(2,2,4)
            hold on
            scatter(G+p,mRes(i).gcTS_rel,5,'r');
            xlabel('G')
            ylabel('STS')
            % Variance
            mRes(i).gcVar = sum(nanvar(gc'));
            mRes(i).mfVar = sum(nanvar(input'));
            i = i + 1;
        end
        save("Temporal_STS2/rX_STS" + ntr+ "trials_3000popgc_50mf_typeOU_ouTarget_t1000_gaus" + gw + "_pop" + p + "_seed_" + randi([1 1000000],1,1) + "_seed2_" + randi([1 1000000],1,1) + ".mat")
    end
end

