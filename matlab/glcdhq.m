% Direct DM Approximation of Anisotropic Gaussians

% Reference: Hanebeck, Huber, Klumpp, "Dirac Mixture Approximation of Multivariate Gaussian Densities", CDC 2009, 
% with some changes of formula in Theorem III.2 according to Matlab code from Uwe Hanebeck:
% svn checkout svn+ssh://i81server.iar.kit.edu/SVN/Publ/2009/DM_Gauss_Approx 

% Compute DM Approximation of Anisotropic Gaussian
% - x    : (L x 1) column vector of x values, initial guess
% - y    : (L x 1) column vector of y values, initial guess
% - wxy  : (L x 1) sample weights
% - SX   : (1 x 1) x standard dev
% - SY   : (1 x 1) y standard dev
%
% Example:
% >> L=40; [x,y]=glcdhq( rand(L,1), rand(L,1), ones(L,1)/L, 1, 1); figure; scatter(x,y); axis equal 
%
function [x, y, G] = glcdhq(x,y, wxy, SX, SY)

    assert(iscolumn(x), 'x must be column vector')
    assert(iscolumn(y), 'y must be column vector')
    assert(iscolumn(wxy), 'wxy_ must be column vector')

    G1 = 0;
    flag = false;
    bmax = 100;
    
    xb=Marshal(x, y);
    %options=optimset('GradObj','on', 'Display','on', 'TolX',1e-22, 'TolFun',1e-22, 'MaxFunEvals',1000000, 'MaxIter',5000, 'LargeScale','off', 'Algorithm','quasi-newton', 'HessUpdate','lbfgs', 'DerivativeCheck','off'); %
    options=optimoptions('fminunc', 'SpecifyObjectiveGradient',true, 'Display','off', 'StepTolerance',1e-22, 'OptimalityTolerance',1e-22, 'ObjectiveLimit',1e-22, 'MaxFunctionEvaluations',1000000, 'MaxIterations',5000, 'Algorithm','quasi-newton', 'HessUpdate','bfgs', 'CheckGradients',false); %
    FunctionWrapper(xb); % test
    [xb, G]=fminunc(@FunctionWrapper, xb, options);
    [x, y]=UnMarshal(xb);


    % Wrapper for distance measure and gradient
    function [f, r]=FunctionWrapper(xb)
        [x, y]=UnMarshal(xb);
        lambda=[1000; 1000; 10; 10; 10]; % moment matching weights 
        f=Distance_Measure_Gaussian_numeric(wxy, x, y, SX, SY, lambda);
        r=Gradient_Gaussian_numeric(wxy, x, y, SX, SY, lambda);
    end


    % Distance measure for axis-aligned Gaussian with different variances
    function G=Distance_Measure_Gaussian_numeric(wxy, x, y, SX, SY, lambda)
        sx=SX;
        sy=SY;

        Cb=log(4*bmax^2)-0.577216;
        b=linspace(0.0001, sqrt(bmax), 100).^2;

        xx=[x,y]; % emulation of new data structure
        N=2;

        L=length(wxy);

        % First Part
        if ~flag
            G1=pi^(N/2)*b.^(N+1)./( sqrt(sx^2+b.^2).*sqrt(sy^2+b.^2) );
            pp=spline(b, G1);
            G1=integral(@(x)ppval(pp,x), 0, bmax);
            flag=true;
        end

        % Second Part
        G2t=0;

        for i=1:L
            G2t=G2t+wxy(i)*exp(-0.5*( xx(i,1)^2./(sx^2+2*b.^2) + xx(i,2)^2./(sy^2+2*b.^2) ));
        end

        G2=(2*pi)^(N/2)*b.^(N+1)./( sqrt(sx^2+2*b.^2).*sqrt(sy^2+2*b.^2) ) .* G2t;

        pp=spline(b, G2);
        G2=integral(@(x)ppval(pp,x), 0, bmax);

        % Alternative Third Part
        Mxx=bsxfun(@minus, x', x);
        Myy=bsxfun(@minus, y', y);
        T=Mxx.^2+Myy.^2;
        %G3=sum(sum( (wxy*wxy') .* xlog(T)));
        G3=pi/8 * wxy'*( 4*bmax^2*exp(-0.5*T/(2*bmax^2)) - Cb*T + xlog(T) -T.^2/(4*bmax^2) )*wxy;
        
        % Combine G1, G2, G3
        G = G1 -2*G2 + G3; 

        % Moment Matching (see also in Gradient_Gaussian_numeric)
        %G = G + lambda(1)*(wxy'*x)^2 + lambda(2)*(wxy'*y)^2 + lambda(3)*(wxy'*(x.^2)-sx^2)^2 + lambda(4)*( wxy'*(x.*y) )^2 + lambda(5)*(wxy'*(y.^2)-sy^2)^2;
    end

    % Gradient for axis-aligned Gaussian with different variances
    function Gr=Gradient_Gaussian_numeric(wxy, x, y, SX, SY, lambda)
        s(1)=SX;
        s(2)=SY;

        Cb=log(4*bmax^2)-0.577216;

        xx=[x,y]; % emulation of new data structure
        N=2;

        L=length(wxy);

        % First part
        db=0.005;
        b=db:db:bmax;

        H=2*(2*pi)^(N/2)*b.^(N+1)./( sqrt(s(1)^2+2*b.^2) .* sqrt(s(2)^2+2*b.^2) );

        Gr1 = NaN(2*L,length(b));
        for eta=1:2
            k=H./( (s(eta)^2+2*b.^2) );
            for i=1:L
                Gr1((eta-1)*L+i, :) = wxy(i) * xx(i, eta) * k .* exp(-0.5*( xx(i,1)^2./(s(1)^2+2*b.^2) + xx(i,2)^2./(s(2)^2+2*b.^2) ));
            end
        end

        Gr1=db*sum(Gr1,2);

        % Second part
        Mxx=bsxfun(@minus, x', x);
        Myy=bsxfun(@minus, y', y);

        M=Mxx.^2+Myy.^2;
        T=plog(M) - M/(4*bmax^2);

        rx=(wxy'*( Mxx.*T ))';
        ry=(wxy'*( Myy.*T ))';

        Gr2=pi*[rx + Cb*(wxy'*x-x); ry + Cb*(wxy'*y-y)]/(2*L);

        Gr3=[2*lambda(1)*(wxy'*x)*wxy + 4*(wxy'*(x.^2)-s(1)^2)*lambda(3)*wxy.*x + 2*( wxy'*(x.*y) )*lambda(4)*wxy.*y; ...
            2*lambda(2)*(wxy'*y)*wxy + 4*(wxy'*(y.^2)-s(2)^2)*lambda(5)*wxy.*y + 2*( wxy'*(x.*y) )*lambda(4)*wxy.*x];

        Gr=Gr1+Gr2;

        % Moment Matching (see also in Distance_Measure_Gaussian_numeric)
        % Gr = Gr + Gr3; 
    end


% Various log-functions
    function y = plog(x)
        indx = (x==0);
        zero_vals = x(indx);
        x(indx) = ones(size(zero_vals));
        y = reallog(x);
    end


    function [x, y] = UnMarshal(xb)
        L = (length(xb))/2;
        x = xb(1:L);
        y = xb(L+1:2*L);
    end


    function xb = Marshal(x, y)
        xb = [x; y];
    end


    function y = xlog(x)
        indx = (x==0);
        zero_vals = x(indx);
        x(indx) = ones(size(zero_vals));
        y = x.*reallog(x);
    end


end