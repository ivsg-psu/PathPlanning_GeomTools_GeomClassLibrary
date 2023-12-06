close all

while(1)

    for ith_frame = 1:48
        drawframe(ith_frame);
        pause(0.01);
    end
end

% https://www.mathworks.com/matlabcentral/communitycontests/contests/6/entries/13907


function drawframe(frame)
    % Get the current global time in UTC
    currentTime=datetime('now');
    
    % Extract hours, minutes, and seconds from the current time
    hours=hour(currentTime);
    hours=mod(hours,12);
    minutes=minute(currentTime);
    seconds=second(currentTime);
     
    % Draw clock design
    cla
    a=[4:12,1:3];
    for i=1:12
        polarplot([pi/6,pi/6]*i,[1.8,2],'LineWidth',4,'Color','k')
        hold on
        text(pi/6*i,1.5,num2str(a(i)),'FontSize',20,'HorizontalAlignment','center')
    end    
    for i=1:60
        polarplot([pi/30,pi/30]*i,[1.9 2],'LineWidth',2,'Color','k')
        hold on
    end
    
    % Calculate angles for the clock hands
    hourAngle=deg2rad(30*(hours+minutes/60)-90);
    minuteAngle=deg2rad(6*(minutes+seconds/60)-90);
    secondAngle=deg2rad(6*seconds-90);
    
    % Draw hour hand
    polarplot([0,hourAngle],[0,1],'LineWidth',6,'Color','k');
    
    % Draw minute hand
    polarplot([0,minuteAngle],[0,1.3],'LineWidth',4,'Color','k');
    
    % Draw second hand
    polarplot([0,secondAngle],[0,1.7],'LineWidth',2,'Color','r');
    
    % Set axis properties
    ax=gca;
    ax.RTickLabel={};
    ax.ThetaTickLabel={};
    ax.ThetaDir='clockwise';
    
    % Set the background color of the figure to black
    l=gcf;
    l.Color='k';
    
    % Turn off grid
    grid off
    
    % Display current time with timezone as title
    title(char(currentTime)+" UTC",'Color','w','FontSize',13);
end