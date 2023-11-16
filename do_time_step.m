function [v, i_impr] = do_time_step(t, v, vb, k, sigma_t, tau_rec, tau_s)
%DO_TIME_STEP Summary of this function goes here
%   Detailed explanation goes here
    arguments
        t       (1,:) double {mustBeReal}
        v       (1,:) double {mustBeReal}
        vb      (1,1) double {mustBeReal}
        k       (1,1) double {mustBeReal}
        sigma_t (1,1) double {mustBeReal}
        tau_rec (1,1) double {mustBeReal}
        tau_s   (1,1) double {mustBeReal}
    end

    m = find(isnan(v), 1, 'first');
    dt = t(2) - t(1);

    fm_impr = NaN(1, m);
    for n = 1 : 1 : m
        fm_impr(n) = sum( exp(- (t(m) - t(n : m)) / tau_s) );
%         hm(n) = exp(- 0.5 * (t(n) / sigma_t) .^ 2) .* exp(- (t(m) - t(n)) / tau_rec) * k * dt ^ 2;

%         fm_impr = NaN(1, m - n + 1);
%         for q = n : 1 : m
%             fm_impr(q) = exp(- (t(m) - t(q)) / tau_s);
%         end

%         Fm_impr(n) = hm(n) * sum(fm_impr);
    end
%     hm = exp(- 0.5 * (t(1 : m) / sigma_t) .^ 2) .* exp(- (t(m) - t(1 : m)) / tau_rec) * k * dt ^ 2;

    F_impr = sum( (k * dt ^ 2) * exp(- 0.5 * (t(1 : m) / sigma_t) .^ 2) ...
        .* exp(- (t(m) - t(1 : m)) / tau_rec) .* fm_impr );
    i_impr = vb * F_impr;

    v(m) = 0;
end

