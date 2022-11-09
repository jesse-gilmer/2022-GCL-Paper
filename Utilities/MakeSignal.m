function [output,stats] = MakeSignal(samples,channels,type,options)
%% Definitions:
%% Samples:
% The number of timepoints generated.

%% :
% The number of input , or, MFs.


%% Types:
% 1 = OU Process
% 2 = Cosine
% 3 = Pulse
% 4 = Step


%% Get output from type.
stats = struct();

switch(type)
    case 1
        mu = options.mu;
        sigma = options.sigma;
        tau = options.tau;
        output = MakeOU(samples,channels,mu,sigma,tau);

        stats.type = 'OU';
        
    case 2
        d = dir("ProcessedDatasets");
        M = load("ProcessedDatasets/"+d(randi([3 length(d)],1,1)).name);
        w = gausswin(25);
        w = w/sum(w);
        Mus = filtfilt(w,1,M.Muscle')';
%         for i = 1:size(Mus,1)
%             Mus(i,:) = Mus(i,:) ./ std(Mus(:)) .* mean(Mus(:)) ;
% %             Mus(i,:) = (Mus(i,:) - min(Mus(i,:)))/(max(Mus(i,:)) - min(Mus(i,:)));
%         end
        Mus(isnan(Mus)) = 0;
%         Mus = (Mus' - min(Mus'))./max(Mus') - min(Mus');
        Mus = Mus' ./ nanstd(Mus') ;
        Mus = Mus ./ nanmean(Mus);
        M0 = mean(Mus');
        Mus(M0<=0,:) = [];
        output = Mus';
        
        
    case 3
        
        
    case 4
        
end

%% Write down universal stats:
stats.numerType = type;
stats.samples = samples;
stats.channels = channels;
     
%% Maker functions:
function [signal] = MakeOU(samples,channels,mu,sigma,tau)
    dt = 1;
    nPts = samples;
    r0 = randn(channels,1)*sigma;
    mf = zeros(channels, nPts) + r0;
    
    for i = 2:nPts
        mf(:,i) = mf(:,i-1)*exp(-dt/tau) + sigma*sqrt(1-exp(-2*dt/tau))*randn(channels,1);
    end
    mf = mf + mu;
    signal = mf;
end

end
