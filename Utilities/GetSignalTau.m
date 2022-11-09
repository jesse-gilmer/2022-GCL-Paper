function [tau] = GetSignalTau(signal)

E = sum(signal);
try
signal(:,E==0) = [];
end

dt = 1e-3;
t = size(signal,1);
lags = 0:dt:(size(signal,1)/1000)-dt;
acs = [];
for n = 1:size(signal,2)
    [ac s] = xcov(signal(:,n));
    acs(n,:) = ac(s>=0)./length(t);
end

% opts = fitoptions('exp1');
P = nanmean(acs);
P = P - min(P);
% % P(P<0) = 0;
% f = fit( lags', P', 'exp1','StartPoint',[1 1],'upper',[Inf 0]);
% tau = round(-100/f.b)/100;
% % plot(f)
% % hold on
% % plot(lags,P)

X  = lags;
Y = P;
fontSize = 14;
% Convert X and Y into a table, which is the form fitnlm() likes the input data to be in.
tbl = table(X', Y');
% Define the model as Y = a + exp(-b*x)
% Note how this "x" of modelfun is related to big X and big Y.
% x((:, 1) is actually X and x(:, 2) is actually Y - the first and second columns of the table.
modelfun = @(b,x) b(1) * exp(x(:, 1)*-b(2));  
beta0 = [P(1), ((P(1)-P(2))/2)]; % Guess values to start with.  Just make your best guess.
% Now the next line is where the actual model computation is done.
mdl = fitnlm(tbl, modelfun, beta0);
% Now the model creation is done and the coefficients have been determined.
% YAY!!!!

% Extract the coefficient values from the the model object.
% The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
coefficients = mdl.Coefficients{:, 'Estimate'};
% Create smoothed/regressed data using the model:
yFitted = coefficients(1) * exp(-coefficients(2)*X);
% % Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
% figure(1)
% subplot(1,2,2);
% hold on;
% plot(X, yFitted, 'r-', 'LineWidth', 2);
% scatter(X,Y);
% grid on;
% title('Exponential Regression with fitnlm()', 'FontSize', fontSize);
% xlabel('X', 'FontSize', fontSize);
% ylabel('Y', 'FontSize', fontSize);
% legendHandle = legend('Noisy Y', 'Fitted Y', 'Location', 'north');
% legendHandle.FontSize = 30;
% formulaString = sprintf('Y = %.3f * exp(-%.3f * X)', coefficients(1), coefficients(2))
% text(7, 11, formulaString, 'FontSize', 25, 'FontWeight', 'bold');
% 
% % Set up figure properties:
% % Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % Get rid of tool bar and pulldown menus that are along top of figure.
% % set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% % Give a name to the title bar.
% set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 

tau = (1000/coefficients(2));

end


