function [J,Icorrect] = CorrectInhomogeneity(I,radius,varargin)
% construct rough estimate of the signal inhomogeneity via morphological closing

% handle variable inputs
if nargin == 3
    filter = varargin{1};
else
    filter = [];
end

% estimate inhomogeneity field
Iinhomogeneity = imclose(I,strel('disk',radius));

% imshow(Iinhomogeneity,[])

% construct the correction field by the estimated signal inhomogeneity
Icorrection = (mean(mean(Iinhomogeneity)))./double(Iinhomogeneity);
Icorrect = double(I).*Icorrection;
if ~isempty(filter)
    if strcmp(filter,'medfilt')
        Icorrect = medfilt2(Icorrect,[3,3]);
    else
        warning('unsupported option for filter provided')
    end
end
% rescale image to the maximum intensity
Imax = max(max(Icorrect));
J = uint8((Icorrect/Imax)*255);
end

