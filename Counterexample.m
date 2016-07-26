% Counterexample for connectedness property of Grady's paper:
% Leo Grady, "Random Walks for Image Segmentation", IEEE Trans. on Pattern 
% Analysis and Machine Intelligence, Vol. 28, No. 11, pp. 1768-1783, 
% Nov., 2006.
%
% 2009.3/4 - Ming-Ming Cheng
% Based on the paper:
% Ming-Ming Cheng, Guo-Xin Zhang, "Connectedness of Random Walk Segmentation",  
% IEEE Trans. on Pattern Analysis and Machine Intelligence, Vol. **, No.
% **, pp. ****-****, ***.,2010
%
% Available at: http://cg.cs.tsinghua.edu.cn/people/~cmm/
%
% Note: Requires installation of the Graph Analysis Toolbox available at:
% http://eslab.bu.edu/software/graphanalysis/
% 

clear;
close all

%% Set experiment parameters
imgS = 512; 
img = 0.9 * ones(imgS, imgS, 3);
inR = 80; outR = 160;
colors = [0.9, 1, 1; 1, 0.9, 1; 1, 1, 0.9;]*0.9;


nTheta = 2;  % number of seeds
r1 = 120; r2 = 240; % seeds positions
%% construct image
for y = 1 : imgS
    for x = 1 : imgS
        tx = x - imgS/2;
        ty = y - imgS/2;
        r = sqrt(tx*tx + ty*ty);
        if r < inR
            img(y, x, :) = [0, 0, 0];
        elseif r > outR
            img(y, x, :) = [0.9, 0.9, 0.9];
        else
            theta = atan2(ty, tx);
            n = (theta / pi + 1)* nTheta ;
            n = floor(n / 2);
            if mod(nTheta, 2) == 0
                img(y, x, :) = colors(mod(n, 2) + 1, :);
            elseif nTheta == 7
                if n == 6
                    img(y, x, :) = colors(3, :);
                else
                    img(y, x, :) = colors(mod(n, 2) + 1, :);
                end
            elseif nTheta < 13
                img(y, x, :) = colors(mod(n, 3) + 1, :);
            else
                disp('Invalidate nTheta not defined currently');
            end
        end
    end
end

imshow(img);

%% calculate seeds and lables
theta = linspace(0, 2 * pi, nTheta + 1); theta = theta(1:end-1) + pi + pi / nTheta;
seeds = zeros(nTheta*2, 2);
seeds(1:nTheta, 1) = r1 * cos(theta) + imgS/2;
seeds(1:nTheta, 2) = r1 * sin(theta ) + imgS/2;
seeds(nTheta + 1:2*nTheta, 1) = r2 * cos(theta) + imgS/2;
seeds(nTheta + 1:2*nTheta, 2) = r2 * sin(theta) + imgS/2;
seeds = ceil(seeds);

labels = [1:2*nTheta];
labels(nTheta + 1 : 2*nTheta) = nTheta+ 1;

%% calculate seed index 
seedsInd = zeros(length(seeds), 1);
for i = 1 : length(seeds)
    seedsInd(i) = sub2ind([imgS imgS], seeds(i, 2), seeds(i, 1)); 
end

%% Apply the random walker algorithms
% [mask,probabilities] = random_walker(img,[sub2ind([X Y],s1y,s1x), ...
%     sub2ind([X Y],s2y,s2x)],[1,2]);
[mask,probabilities] = random_walker(img,seedsInd,labels);

%% Display results
% oringnal figure and seeds
subplot(221);
imagesc(img);
colormap('gray')
axis equal tight off
hold on
plotOpt = ['r.'; 'r+'; 'bx';  'g.'; 'b.'; 'm.'; 'y.'];
for i = 1 : nTheta
    plot(seeds(i, 1), seeds(i,2), plotOpt(mod(i, length(plotOpt))+1, :));
end
plot(seeds(1+nTheta:end, 1), seeds(1+nTheta:end, 2), 'k*');
title('Image with seeds')

% mask
subplot(222);
imagesc(mask)
axis equal tight off
title('Output mask');


% Two typical probality image
subplot(223);
imagesc(probabilities(:,:,1))
axis equal tight off
title('Probability 1st inside seed')

subplot(224);
imagesc(probabilities(:,:,nTheta+1))
axis equal tight off
title('Probability outside seed')

disp(sprintf('n = %d: %f, %f', nTheta, probabilities(imgS/2,imgS/2,1), probabilities(imgS/2,imgS/2,nTheta+1)));