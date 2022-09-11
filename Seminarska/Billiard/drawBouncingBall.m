function drawBouncingBall(trajPts,tableBorder,ballR)
	% function drawBouncingBall(trajectoryPts,tableBorder,ballR)
	% 
	% Prikaze animacijo, v kateri se kroglica odbija od robov biljardne mize po
	% odbojnem zakonu. 
	%
	% Vhod:
	% trajPts - matrika(2xstTock) zaporednih tock na trajektoriji po kateri 
	% potuje kroglica
	% tableBorder - matrika(2xstOglisc+1) tocke, ki predstavljajo oglisca 
	% biljardne mize
	% ballR - radij kroglice, ki se odbija po biljardni mizi

	% novo okno za diagrame
	f = figure;

	% izklopimo nepotrebne stvari
	f.MenuBar = "none";
	axis equal;
	axis off;

	% meje grafikonov fiksiramo
	xlim([-1.04 1.04]);
	ylim([-1.04 1.04]);

	hold on;
	% narisemo nortanji(od katerega se kroglica odbija) in zunanji rob(za lepsi izgled)
	plot(tableBorder(1,:),tableBorder(2,:),"b");  
	plot(tableBorder(1,:)*1.03,tableBorder(2,:)*1.03,"b");

	circleApproxN = 60;  % 60 tock s katerimi prametriziramo kroznico
	ballParam = linspace(0,2*pi,circleApproxN); 
	ball = [ballR*cos(ballParam);ballR*sin(ballParam)]; % kroznica z radijem ballR
	numFrames = size(trajPts,2); % stevilo slicic, ki jih moramo animirati

	% koliko slicic na cekundo zelimo prikazati
	fps = 60;
	refreshRate = 1/fps;

	% grafikon katerega tocke bomo posodabljali 
	bouncingBall = plot(0:circleApproxN:0,0:circleApproxN:0,"k");
	numTrails = 20;
	
	% sled
	trail = plot(0:numTrails:0,0:numTrails:0,"r:");
	trail.LineWidth = 2.5;

	% v zivo posodabljamo podatke za kroglico, kar daje obcutek gibanja
	for i = 1:numFrames
		% sredisce kroznice premikamo na koordinate trajektornih tock
		bouncingBall.XData = ball(1,:) + trajPts(1,i);
		bouncingBall.YData = ball(2,:) + trajPts(2,i);

		% sled risemo le ko smo ze prepotovali dolzino sledi
		if i >= numTrails
			trail.XData = trajPts(1,(i-numTrails+1):i);
			trail.YData = trajPts(2,(i-numTrails+1):i);
		end

		pause(refreshRate); % zagotovimo fps slicic na sekundo
	end
	hold off;
end