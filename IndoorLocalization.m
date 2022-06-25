
N = 3;  % number of anchors
M = 5;  % number of users

% distance dependent err SNR
distMeasurementErrRatio = 0.1;  % it means that the accuracy of distance measurement is 90 %
								% for instance the inaccuracy of a 1m measured distance
								% is around .1 meter.

networkSize = 50;  % we consider a 50x50 area that the mobile can wander

anchorLoc   = [0  0; % set the anchor at 4 vertices of the region
			   -10  0;
			   0   -10];

% building a random location for the mobile node
%mobileLoc  = networkSize*rand(M,2);

mobileLoc = [15 0;
			20 5;
			25 20;
			27 28;
			35 32];


% Computing the Euclidian distances    
	% very fast computation :)
	% distance   = sqrt(sum( (anchorLoc - repmat(mobileLoc,N,1)).^2 , 2));
	
% easy to understand computation
distance = zeros(N,M);
for m = 1 : M
	for n = 1 : N
			distance(n,m) = sqrt( (anchorLoc(n,1)-mobileLoc(m,1)).^2 + ...
										(anchorLoc(n,2)-mobileLoc(m,2)).^2  );
	end
end
% Plot the scenario
f1 = figure(1);
clf
% plot(anchorLoc(:,1),anchorLoc(:,2),'ko','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
grid on
hold on
plot(mobileLoc(:,1),mobileLoc(:,2),'b+','MarkerSize',8,'lineWidth',2);

% noisy measurements
distanceNoisy = distance + distance.*distMeasurementErrRatio.*(rand(N,M)-1/2);

% using gaussian 


numOfIteration = 5;

% Initial guess (random locatio)
mobileLocEst = networkSize*rand(M,2);
% repeatation
for m = 1 : M
	for i = 1 : numOfIteration
		% computing the esimated distances
		distanceEst   = sqrt(sum( (anchorLoc - repmat(mobileLocEst(m,:),N,1)).^2 , 2));
		% computing the derivatives
			% d0 = sqrt( (x-x0)^2 + (y-y0)^2 )
			% derivatives -> d(d0)/dx = (x-x0)/d0
			% derivatives -> d(d0)/dy = (y-y0)/d0
		distanceDrv   = [(mobileLocEst(m,1)-anchorLoc(:,1))./distanceEst ... % x-coordinate
						 (mobileLocEst(m,2)-anchorLoc(:,2))./distanceEst];   % y-coordinate
		% delta 
		delta = - (distanceDrv.'*distanceDrv)^-1*distanceDrv.' * (distanceEst - distanceNoisy(:,m));
		% Updating the estimation
		mobileLocEst(m,:) = mobileLocEst(m,:) + delta.';
	end
end    
plot(mobileLocEst(:,1),mobileLocEst(:,2),'ro','MarkerSize',8,'lineWidth',2);
legend('True location','Estimated location',...
	   'Location','Best')

% Compute the Root Mean Squred Error
Err = mean(sqrt(sum((mobileLocEst-mobileLoc).^2)));
xlabel('X Position (m)');
ylabel('Y Position (m)');
ax = gca;
ax.FontSize = 16;
title({'Tx power: 23dBm;', strcat('Mean Estimation Error: ',num2str(Err),'meters')})
axis([-0.1 1.1 -0.1 1.1]*networkSize)