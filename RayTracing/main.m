% Ray Tracing for Indoor Surveillance
% Inheritted from the code of Salaheddin Hosseinzadeh
%       1- Free Space Loss under 1 meter (done)
%       2- Reflection and refraction coefficients (done)
%       3- Antenna Gain (done)
%       4- Angle of incidence (done)
%       5- Second reflections (done)
%       6- Structural design import from AutoCadv(done)
%       7- Polarization (done)
%       8- Angle of incidence (done)
%       9- Angle of incidence and arrival (done)
%       10- Fresnel Coeffs (done)
%

% There are two ways to define walls or structure
% 1- Method one (easy 2D to 3D):
%       This is the easy way of doing it. If the walls are all of the same
%       height (ususally are), ignore their height and think of them as
%       lines to make it 2D. Lets say walls height is 4 metters, then
%       change "ceilingLevel = 4".
%       -So, now walls are lines, and they have 2
%       ends. You can define them as a CSV file now. Each wall has a start
%       and a end defined by [Xstart,Ystart,Zstart] and [Xend,Yend,Zend].
%       -Ignore Zstart and Zend, they are always zero. 
%       -Lets say the reference point of start is [0,0,0], this is where
%       the first wall starts from, and it end point is at [36,0,0]. The
%       first line in the excel or CSV file then should look like this.
%       0  ;0  ;0     
%       36 ;0  ;0 


clear all
clc

%% Ray Tracing Engine Parameters    

optimizationMode = 0;   % if 1, then reduces the Rx mesh size to measurement locations only
plotMode = 1;           % Ray Tracing engine plot mode
demoMode = 0;           % Ray Tracing Engine Plot Mode

losFlag = 1; 
reflectionFlag = 2;                     % whether or not calculate First reflections
secondReflectionFlag = 1;
reflectExaggerationFac = 1e0;          % Must be 1 unless to emphasize reflection for demonstration purposes
                                    % whether or not calculate LoS
disableIncidentAngle = 0;               % 1 Disables the incident angle calculation, if disableIncidentAngle= 1
solidIncidentAngle = 45;                % if disableIncidentAngle =1, then assign this which overwrites all the incident angles! This is unnecessary feature
polarizationSwap = 1;               % (See notes in HOW THIS WORSK)  % 1, Applies TE to walls and TM to ceiling. 0, applies TM to the walls and TE to the ceiling


imageRSSIScale = 5;         % increase this if number of meshes nodes are small
grayScaleImage = 0;


freq = 100e6;  % frequency in hertz
lightVel = 3e8;
lambda = lightVel./freq;
refDistance = 1;            % Reference distance from Tx
FPSLRefLoss = 0;


antennaGainRes = 40;
antennaEffiLoss = -11.5;         % dB antenna efficiency, missmatch and loss all together

ceilingEnable = 0; % Allowing to define ceiling and floor
groundLevel = 0;
ceilingLevel = 4;  % Height of the ceiling

% 
mesh_.xNodeNum = 40;   % Keep the x and y mesh size the same, increase the size for better resolution and especially if you're increasing the frequency
mesh_.yNodeNum = 40;
mesh_.zNodeNum = 1;


%% Antenna Gain pattern calculation

[TxAntennaGainAE] = AntennaTemp (antennaGainRes,demoMode) + antennaEffiLoss;  % TxAntennaGainAE needs to be in dB
RxAntennaGainAE = TxAntennaGainAE;


% Location of the transmitter (s)

Tx.xyz = [
            35,32 ,0
%            18,10,1.5
                ];
            
% power of the transmitter dB(m)
Tx.power =  [
             -23
%            -10
                ]; % Ray Power at 1m in dB

% Defining the boundary of the analysis (something like a boundary condition) 
boundary = [
            -5,50
            -3,50
            -0,3
            ];    



% Reads the structure from an excel file (see in this code section at the
% top)
[wallxyz1, wallxyz2, wallxyz3, wallxyz4,wallX,wallY,wallZ] = CSV23D_V1(demoMode,groundLevel,ceilingLevel,Tx.xyz);

wall.xyz1 = wallxyz1;
wall.xyz2 = wallxyz2;
wall.xyz3 = wallxyz3;
wall.xyz4 = wallxyz4;

wall.X = wallX;
wall.Y = wallY;
wall.Z = wallZ;


% Define the ceiling of the structure manually if required walls can be
% defined the same fashion
if ceilingEnable == 1

    ceillFloor.xyz1 = [0,0,ceilingLevel
                       50,35,ceilingLevel
                       0,0,groundLevel
                        ];

    ceillFloor.xyz2 = [0,50,ceilingLevel
                        50,35,ceilingLevel
                        0,30,groundLevel

                        ];

    ceillFloor.xyz3 = [50,23.9,ceilingLevel
                        50,27.83,ceilingLevel
                        50,27.54,groundLevel
                        ];

    ceillFloor.xyz4 = [50,0,ceilingLevel
                        50,23.9,ceilingLevel
                        50,0,groundLevel
                        ];
else
    
    ceillFloor.xyz1 = [];
    ceillFloor.xyz2 = [];
    ceillFloor.xyz3 = [];
    ceillFloor.xyz4 = [];
                    
end

% Relative permittivity of falls can be defined here individually
wall.relativePerm = 6*ones(size(wall.xyz1,1)+size(ceillFloor.xyz1,1),1);



%% Adding Ceillilng and Floor to the structure
for i = 1:size(ceillFloor.xyz1,1)
    
    wall.xyz1 = [wall.xyz1;ceillFloor.xyz1(i,:)];
    wall.xyz2 = [wall.xyz2;ceillFloor.xyz2(i,:)];
    wall.xyz3 = [wall.xyz3;ceillFloor.xyz3(i,:)];
    wall.xyz4 = [wall.xyz4;ceillFloor.xyz4(i,:)];
    
end

RayLib



