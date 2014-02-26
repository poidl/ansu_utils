% Load some image data.
load clown
% Define a colormap that consists of two separate colormaps.
% The object in the left axes will use the first half of the
% colormap, and the object in the right axes will use the
% second half.
cmap1 = map;
cmap2 = gray(size(map,1));
cmap = [cmap1;cmap2];
colormap(cmap)
% Generate the first image.
subplot(121)
image(X)
% Generate the second image.  In this axes, the image data 
% is scaled to start at 1+<the maximum value of other 
% image>.  This is the equivalent of setting the CData 
% for each image so that they have contiguous, nonoverlapping 
% values.
subplot(122)
X2 = X + max(X(:));
image(X2)
%Note that if surface or patch plots were used instead of images, you would also have to set the CLim properties of each axis to [min(X(:)) max(X2(:))]. The following code fragment illustrates how to do this using the SET function rather than CAXIS:

ax = findobj(gcf,'Type','axes');
set(ax,'CLim', [min(X(:)) max(X2(:))])
%Note that we can set the CLim property of both axes simultaneously with SET; CAXIS will only affect one axis.

%If the CData is logical, it should be converted to direct indexing in order to use more than the first 2 rows of the colormap. For example, following the same process as above, but with a binary image:

load clown;
%Create logical CData 
X = X>mean(X(:));
subplot(121)
h(1) = image(X);
subplot(122)
h(2) = image(X);
colorbar
cmap1 = summer(2);
cmap2 = copper(2);
cmap = [cmap1;cmap2];
colormap(cmap)
%CData for h(1) and h(2), should be converted from logical 
%indexing to direct indexing in order to make use of more than 2 colors  
%"1" is of type double by default, so adding to logical CData causes 
%the result to be converted to double due to class precedence
h1CData = get(h(1),'CData');
%This converts from logical to direct indexing
h1CData = h1CData + 1;
set(h(1),'CData', h1CData );
h2CData = get(h(2),'CData');
%First convert from logical to direct indexing
%Conversion from logical can also be done directly using DOUBLE:
h2CData = double(h2CData);
%Convert to direct index starting at 1 instead of 0
h2CData = h2CData + 1;
%Then set to different part of colormap than first image
h2CData = h2CData + max(h1CData(:));
set(h(2),'CData', h2CData );