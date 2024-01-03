%% LCD Gaussian Sampling
%  High quality (and not just standard normal reference samples that are transformed) 


%% Dependencies
% [flist,plist] = matlab.codetools.requiredFilesAndProducts('main.m'); [flist'; {plist.Name}']
% glcdhq.m; Optimization Toolbox


%% Sampling

% Gaussian Parameters
L = 30; 
sx = .5;
sy = 1;

% Initialization 
x0 = rand(L,1)*sx;
y0 = rand(L,1)*sy;
w = ones(L,1)/L;

% Compute
[x,y] = glcdhq( rand(L,1), rand(L,1), w, sx, sy); 


%% Plot
fig = figure(292349343);
hs = scatter(x,y); 
set(hs, 'Marker','.', 'SizeData',1000)
axis equal 


exportgraphics(fig, 'gaussian.pdf')
