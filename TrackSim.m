% Track Simulator
% 8/24/17
% Chris Carter

format shortg
% FTireTransMax = MuTire*FTireVert;
% TireSlipRatio = (TireAngVel * TireRadius / Velocity - 1) * 100;

% Sections are [Arc or Straight Distance, Corner Inside Radius] in meters
Track = [1 1000 0 ; 2 15 20]; 
Direction(1) = 0; % Degrees
Forcex = 50; % N
Forcey = 2; % N
mass = 250; % kg
Distancex(1) = 0; % m
Distancey(1) = 0; % m
Velocityx(1) = 0; % m/s
Velocityy(1) = 0; % m/s
Accelerationx = Forcex/mass; % m/s
Accelerationy = Forcey/mass; % m/s
time(1) = 0; %s
DeltaTime = 0.01; % seconds
SectionsDistance = 0;
TotalDistance = sum(Track(:,2));
CurrentDistance = 0;
i = 1;
j = 1;
[Sections,Dummy] = size(Track);

while CurrentDistance < TotalDistance
    i = i + 1;
    time(i) = time(i-1) + DeltaTime;
    Velocityx(i) = Velocityx(i-1) + Accelerationx*(time(i)-time(i-1));
    Distancex(i) = Distancex(i-1) + Velocityx(i) * (time(i)-time(i-1)) + 1/2 * Accelerationx * (time(i)-time(i-1))^2;
    Velocityy(i) = Velocityy(i-1) + Accelerationy*(time(i)-time(i-1));
    Distancey(i) = Distancey(i-1) + Velocityy(i) * (time(i)-time(i-1)) + 1/2 * Accelerationy * (time(i)-time(i-1))^2;
    CurrentDistance = sqrt((Distancex(i) - Distancex(i-1))^2 + (Distancey(i) - Distancey(i-1))^2) + CurrentDistance;
    % Angle in degrees at the current section of track
    if Track(j,3) == 0
        Direction(i) = Direction(i-1);
    else
        Direction(i) = ((CurrentDistance - SectionsDistance) / Track(j,3))*180/ pi;
    end
    % If the current section of track has been covered, move on to the next
    % section of track
    if Track(j,2) <= CurrentDistance - SectionsDistance
        if j == Sections
            break;
        end
        SectionsDistance = SectionsDistance + Track(j,2);
        j = j + 1;
    end
end

Velocity = sqrt(Velocityx(i)^2 + Velocityy(i)^2)
Distance = sqrt(Distancex(i)^2 + Distancey(i)^2)
CurrentDistance
TotalDistance
Accuracy = abs((CurrentDistance-TotalDistance)/TotalDistance)*100
Direction = Direction(i)
