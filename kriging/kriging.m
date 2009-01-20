function [z, varargout] = kriging(xi, yi, zi, x, y, pl)
% Perform interpolation using Ordinary kriging.
% [z, z_var] = kriging(xi, yi, zi, x, y) 
%
% [x, y, z] = kriging(...)
% [x, y, z, z_var] = kriging(...)
% [x, y, z, z_var, time] = kriging(...)
%
% Matt Foster <ee1mpf@bath.ac.uk>

% Code heavily based on recipe from MATLAB recipes for Earth sciences. Martin
% H. Trauth, Springer 2006
    
    warning off all;

    tic;

    if nargin < 6
        pl = 0;
    end

    % Create Variogram -- warning.. even this is slow!

    [X1,X2] = meshgrid(xi);
    [Y1,Y2] = meshgrid(yi);
    [Z1,Z2] = meshgrid(zi);

    D = sqrt((X1 - X2).^2 + (Y1 - Y2).^2);

    G = 0.5*(Z1 - Z2).^2; 

    D2 = D.*(diag(xi*NaN)+1);
    lag = mean(min(D2));
    hmd = max(D(:))/2;
    max_lags = floor(hmd/lag);
    LAGS = ceil(D/lag);

    for i = 1:max_lags
            SEL = (LAGS == i);
            DE(i) = mean(mean(D(SEL)));
            PN(i) = sum(sum(SEL == 1))/2;
            GE(i) = mean(mean(G(SEL)));
    end
    lags=0:max(DE);


    mod_func = fit_spherical(DE, GE, var(zi(:)));

    model = mod_func(D);

    if pl
        plot(DE, GE, '.');
        hold on
        plot(DE, mod_func(DE));
        hold off

        dlmwrite('variogram.dat', [DE', GE', mod_func(DE)'], '\t');

    end

    % add zeros as the lat for and col of the variogram model.
    n = length(xi);
    model(:,n+1) = 1;
    model(n+1,:) = 1;
    % add a zero to the bottom right corner
    model(n+1,n+1) = 0;

    % Invert the variogram model
    model_inv = inv(model);

    Zg = nan(size(x));
    s2_k = nan(size(x));

    % Perform the Kriging
    for k = 1:length(Zg(:));
        DOR = ((xi - x(k)).^2+(yi - y(k)).^2).^0.5;
        G_R = mod_func(DOR);
        G_R(n+1) = 1; 
        E = model_inv*G_R; 
        Zg(k) = sum(E(1:n,1).*zi); 
        s2_k(k) = sum(E(1:n,1).*G_R(1:n,1))+E(n+1,1); 
    end

    whos z
    z = Zg;
    z_var = s2_k;

    time = toc;

    if nargout == 2 || nargout == 3
      varargout{1} = z_var;
    end

    if nargout == 3
      varargout{2} = time;
    end

    if nargout >= 4
      tmp = z;
      z = xx;
      varargout{1} = yy;
      varargout{2} = tmp;
      varargout{3} = z_var;
    end

    if nargout == 5
      varargout{4} = time;
    end
 
    warning on all;

end

function mod_func = fit_spherical(DE, GE, var_z)

    % Only fit to the section below 90% of the variance
    ind = 1:min(find(GE > var_z * 0.9));
    coef = [DE(ind);DE(ind).^3]' \ GE(ind)';
    c = sqrt((-0.5 .* (1.5 ./ coef(1)).^-3)./coef(2));
    a = 1.5 .* c / coef(1);

%----------------------------------------------------
%     lag_var = DE(min(find(GE > var_z)));
%     if a > lag_var;
%         a = lag_var;
%     end
% 
%     if c > var_z
%         c = mean(GE(min(find(GE > var_z)):end));
%     end
%----------------------------------------------------

    if isreal(a)
        % Return an anonymous function that models the variogram
        mod_func = @(h) (coef(1) .* h + coef(2) .* h.^3) .* (h <= a) ...
        + c.*ones(size(h)) .* (h > a);
    else
        mod_func = @(h) (coef(1) .* h + coef(2) .* h.^3);
    end

end

