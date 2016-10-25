function CenterOfMass = center_of_mass(data)
% Find the center of mass of a 3d or 2d data array.

   % Dimensions of the center of the data
   size_data = size(data);
   
   % Number of dimensions in the array data
   dims = length(size_data);
   
   % Sums the data in all dimensions
   mass = sum(data(:));
   
   % Creates empty center of mass vector with zeros in each dimension
   CenterOfMass = zeros(dims,1);
   
   % Loop through each dimension to find the center of mass for that
   % respective axis
   for dim = 1:dims
       % size_dim is ones except it replaces the current dimension with the
       % size, so [X 1 1], [1 X 1], [1 1 X] would be the three values of it
       % for each dim iteration for 3 dimensions.
       size_dim = ones(1,dims);
       size_dim(dim) = size_data(dim);
       
       % The index matrix will be replicated in all dimensions besides the
       % current dim iteration. Set the current dim to 1 because it will
       % not be replicated.
       replicate = size_data;
       replicate(dim) = 1;
       
       % Creates a matrix the same size as data whose values contain the
       % index of each value in data. Index is with respect to the current
       % dim iteration.
       indices = repmat(reshape(1:size_data(dim),size_dim),replicate);
       
       % Center of mass in the given dimension is the sum of the indices
       % times the data divided by total mass.
       CenterOfMass(dim) = sum(indices(:).*data(:))./mass;
   end

end

