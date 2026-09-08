function [Pvisc, data] = objective_kelvin_voigt_piola(F, Fold, par)
    I3 = eye(3);
    J = det(F);
    if J <= 0
        error('Negative or zero J encountered. Element inverted.');
    end

    Finv = F \ I3;
    FinvT = Finv.';
    Fdot = (F - Fold) / par.dt;
    L = Fdot * Finv;
    D = 0.5 * (L + L.');
    trD = D(1,1) + D(2,2) + D(3,3);
    devD = D - (trD/3) * I3;
    etaBulk = objective_kelvin_voigt_bulk_viscosity(par);

    Tvisc = 2 * par.etaE * devD + etaBulk * trD * I3;
    Pvisc = J * Tvisc * FinvT;

    if nargout > 1
        data = struct('J', J, 'Finv', Finv, 'FinvT', FinvT, ...
            'Fdot', Fdot, 'L', L, 'Tvisc', Tvisc, 'etaBulk', etaBulk);
    end
end

function etaBulk = objective_kelvin_voigt_bulk_viscosity(par)
    etaBulk = par.etaE;
    if isfield(par, 'Ke') && isfield(par, 'Ge') && par.Ge > 0
        etaBulk = par.etaE * par.Ke / par.Ge;
    end
    if isfield(par, 'etaBulkE') && isfinite(par.etaBulkE) && par.etaBulkE > 0
        etaBulk = par.etaBulkE;
    end
end