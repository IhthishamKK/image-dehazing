image_name = 'cityscape'; 
img_hazy = imread(['images/',image_name,'_input.png']);

fid = fopen(['images/',image_name,'_params.txt'],'r');
[C] = textscan(fid,'%s %f');
fclose(fid);
gamma = C{2}(1);

A = reshape(estimate_airlight(im2double(img_hazy).^(gamma)),1,1,3);

disp(A);

figure('Position',[50,50, size(img_hazy,2)*3 , size(img_hazy,1)]);
subplot(1,3,1); imshow(img_hazy);    title('Hazy input')

% Read in the image
img = imread('images\cityscape_input.png');

% Display the image
imshow(img);

% Use imtool to view the color channels
imtool(img);

% Read in the image as a matrix
img = imread('images\forest_input.png');

% Get the red, green, and blue color channels
red_channel = img(:,:,1);
green_channel = img(:,:,2);
blue_channel = img(:,:,3);

% Print the values of the first pixel in each channel
disp(['Red channel value: ' num2str(red_channel(1,1))]);
disp(['Green channel value: ' num2str(green_channel(1,1))]);
disp(['Blue channel value: ' num2str(blue_channel(1,1))]);



