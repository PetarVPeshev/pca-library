function varargout = create_time_arrays(dt, time_lim, out_type)
%CREATE_TIME_ARRAYS Summary of this function goes here
%   Detailed explanation goes here
    arguments
        dt       (1,1) double {mustBePositive}
        time_lim (1,2) double {mustBeNonNan, mustBeIncreasing}
        out_type (1,:) char {mustBeMember(out_type, {'arrays', 'struct'})} = 'arrays'
    end

    if time_lim(1) > 0
        error('LowerTimeLimit:higherThanZero', 'Lower time limit is higher than 0 seconds.');
    end

    t_sim = time_lim(1) : dt : time_lim(2);     % Simulation time vector
    t_res = (0 : 1 : length(t_sim) - 1) * dt;   % Impulse response time vector

    if strcmp(out_type, 'arrays')
        if isequal(t_sim, t_res)
            varargout{1} = t_sim;
        else
            varargout{1} = t_sim;
            varargout{2} = t_res;
        end
    else
        varargout = struct('t_sim', t_sim, 't_res', t_res, 'dt', dt);
    end
end
