function [tilt_angle, rotation_angle] = Find_table_tilt_angle(Bdata, voxelSize, limits,offset,showFlag,Nshow)
% Finds and returns the table tilt angle from a set of threshold images.

% If showFlag is set to true, this is the interval of axial slices between
% showing the slice and the line of the flat part of subjects posterior
if nargin < 6
    Nshow = 20;
end

% Whether or not to show the slices and the determined line
if nargin < 5
    showFlag = true;
end
if isempty(showFlag)
    showFlag = true;
end

% This offset is used for determining the left/right half of the subject's
% body. Imagine there is a plane parallel to the sagittal plane that is
% offset by this parameter. The portion of the body to the left is the
% 'left' half and the other is 'right' half.
% Not entirely sure how this parameter will affect the slope calculation
if nargin < 4
    offset = 10;
end
if isempty(offset)
    offset = 10;
end

% Dimensions for the body data
dims = size(Bdata);

% After calculating the lowest row for each slice, the limit is the slices
% that are used to create a linear regression fit with the row and to
% compensate for the axis tilt.
if nargin < 3
    limits = [1 dims(3)];
end
if isempty(limits)
    limits = [1 dims(3)];
end

% Compose a row of zeros that is like a axial slice
yC = zeros(1,double(dims(3)));
rotation_angles = yC;

% Loop through each axial slice
for i = 1:dims(3)
    % build L/R masks
    BdataL = 0*Bdata(:,:,i);
    BdataR = BdataL;
    
    % Left half of BDataL set to Bdata. Left of Sagittal plane
    BdataL(:,1:floor(dims(2)/2)-offset) = Bdata(:,1:floor(dims(2)/2)-offset,i);
    
    % Right half of BDataR set to BData. Right of Saggittal plane
    BdataR(:,ceil(dims(2)/2)+offset:dims(2)) = Bdata(:,ceil(dims(2)/2)+offset:dims(2),i);
    BdataL = (BdataL>0.5);    
    BdataR = (BdataR>0.5);
    
    % compute coordinates of lowest point
    
    % Add up the rows of BdataL.
    sumBdataL = sum(BdataL,2);
    yL = dims(1);
    stopFlag = false;
    for j = (dims(1)-1):-1:1
        % Starting from the lowest row, if the row ahead of it has no
        % values > 0.5 and the current row does, then set this to the
        % lowest row for the left side
        if sumBdataL(j+1) == 0 && sumBdataL(j) > 0 && ~stopFlag
            yL = j;
            stopFlag = true;
        end
    end
    
    % Get all nonzero values from the lowest row, this is xLs
    % Add up xLs and divide by length to get the average value, which is
    % the average location of where the nonzero values are located
    xLs = double(find(BdataL(yL,:)));
    xL = sum(xLs)/length(xLs);
    
    % Add up the rows of BdataR.
    sumBdataR = sum(BdataR,2);
    yR = dims(1);
    stopFlag = false;
    for j = (dims(1)-1):-1:1
        % Starting from the lowest row, if the row ahead of it has no
        % values > 0.5 and the current row does, then set this to the
        % lowest row for the right side
        if sumBdataR(j+1) == 0 && sumBdataR(j) > 0 && ~stopFlag
            yR = j;
            stopFlag = true;
        end
    end
    xRs = double(find(BdataR(yR,:)));
    xR = sum(xRs)/length(xRs);
    
    % Slope between lowest coordinate in left & right half of plane. If no
    % tilting is present in this slice, this will be zero.
    m = (yR - yL)/(xR - xL);
    
    % Rotation angle is inverse tangent of y over x.
    rotation_angles(i) = atan2d((yR - yL),(xR - xL));
    xC = (double(dims(2))+1)/2.0;
    
    % Interpolate lowest point in the middle by multiplying slope on point
    % (xL, yL) to get (center, yC).
    yC(i) = m*(xC-xL) + yL;
    
    % Show every Nshow slices if showFlag is set
    if mod(i,Nshow) == 0 && showFlag
        figure, 
        imshow(Bdata(:,:,i),[]);
        line([xL,xC],[yL,yC(i)],...
               'LineWidth',3,...
               'Color',[.2,.5,1]); 
        line([xC,xR],[yC(i),yR],...
               'LineWidth',3,...
               'Color',[.2,.5,1]);  
        drawnow;
    end        
end
slices = 1:dims(3);


if showFlag
figure;
plot(slices,yC);
xlabel('Scan slice number')
ylabel('y-coord of center of table')
title('Y-center vs slice')

figure;
plot(slices,rotation_angles);
xlabel('Scan slice number')
ylabel('Rotation angles (degrees)')
title('Rotation angle vs slice')
end

% Get the slices between limits(1) and limits(2) and the lowest y for the
% slices
z = slices(limits(1):limits(2))';
y = yC(limits(1):limits(2))';

% perform linear regression y = beta(1) + beta(2)*z 
Z = [ones(length(z),1),z];
beta = Z\y; 

% Rise in mm, run in mm
rise = voxelSize(1)*beta(2)*(limits(2)+1-limits(1));
run = (limits(2)+1-limits(1))*voxelSize(3);

tilt_angle = atan2d(rise,run);

% extremum correction
ymax = 0;
if tilt_angle > 0
    % Loop through each of the z slices
    for k = z(1):z(end)
        % If the lowest center row is less than ymax
        if yC(k) < ymax
            % Since the tilt angle is greater than zero, each subsequent k
            % should have a larger ymax, if not, it corrects by setting it
            % equal to ymax.
            y(k-z(1)+1) = ymax;
        else
            % Set ymax equal to this value
            ymax = yC(k);
        end
    end
else
    % Loop through each of z slices start from end going to beginning
    for k = z(end):-1:z(1)
        % If the lowest center row is less than ymax
        if yC(k) < ymax
            % Since the tilt angle is less than zero, each subsequent k
            % should have a larger ymax, if not, it corrects by setting it
            % equal to ymax.
            y(k-z(1)+1) = ymax;
        else
            % Set ymax equal to this value
            ymax = yC(k);
        end
    end    
end

% perform linear regression y = beta(1) + beta(2)*z 
% Redo calculations from above
beta = Z\y; 
yFit = Z*beta;

rise = voxelSize(1)*beta(2)*(limits(2)+1-limits(1));
run = (limits(2)+1-limits(1))*voxelSize(3);

% In mm/mm
tilt_angle = atan2d(rise,run);

figure,
scatter(z,y)
hold on
plot(z,yFit)
daspect([voxelSize(1) voxelSize(3) voxelSize(1)])
xlabel('Scan slice number')
ylabel('Table coordinate')
title('Table coordinate vs slice')
grid on

%size(rotation_angles)

% No general trend for rotation_angle, just take the mean of the
% rotation_angle. In degrees
rotation_angle = nanmean(rotation_angles,2);

end

