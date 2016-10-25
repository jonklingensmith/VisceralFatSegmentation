function Bdata = Make_DIXON_threshold_image(Fdata, Wdata)
% Creates a combined threshold image from fat and water saturated DIXON
% data.

% Get the dimension of the fat data
dims = size(Fdata);

% Get largest intensity value from sagittal and coronal planes from fat
% and water signals. This is an axial vector containing largest intensity
% values, then pick smallest value out of the slices.
maxData = min([min((max(max(Fdata)))),min((max(max(Wdata))))]);

% Anything that is over the maxData threshold is capped at maxData.
% Not entirely sure why maxData is the normalizing factor
Fdata(Fdata>maxData) = maxData;
Wdata(Wdata>maxData) = maxData;
maxData1 = double(maxData);

% Loop through each axial slice
for slice = 1:dims(3)
    % Normalize the fat & water signals by maxData. Bdata is sum of fat and
    % water signals.
    % Each slice of fat & water signals is rotated by 90deg and then
    % flipped along vertical axis.
    % Each slice is transformed as below as an example:
    % 1 2 3     9 8 7
    % 4 5 6 --> 6 5 4
    % 7 8 9     3 2 1
    Bdata(:,:,slice) = double(fliplr(rot90(Fdata(:,:,slice)',2)))/maxData1 ...
                  + double(fliplr(rot90(Wdata(:,:,slice)',2)))/maxData1;          
end

% This should always returns maxData1
maxData2 = min([maxData1,min(max(max(Bdata)))*maxData1]);
% Should be multiplying by one basically
Bdata = Bdata*(maxData1/maxData2);
% Normalize Bdata, anything over 1.0 caps to 1.0, anything under 0.1 goes
% to 0.0. I'm not sure why 0.1 is the magic number.
Bdata(Bdata>1.0) = 1.0;
Bdata(Bdata<0.1) = 0.0;

end

