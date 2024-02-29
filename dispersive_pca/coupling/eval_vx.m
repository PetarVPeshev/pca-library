function varargout = eval_vx(i, args)
%EVAL_VX Summary of this function goes here
%   Detailed explanation goes here
% TODO: x_feed argument is not strictly necessary (calculating distance from feed can be done outside routine)
    arguments
        i           (1,:) double            {mustBeNonNan}
        args.x      (1,:) double            {mustBeNonNan, mustBeReal} = double.empty(1, 0)
        args.t_res  (1,:) double            {mustBeNonNan, mustBeReal} = double.empty(1, 0)
        args.f      (1,:) double            {mustBeNonNan, mustBeReal} = double.empty(1, 0)
        args.slot         SlotInDielectrics                            = SlotInDielectrics.empty
        args.x_feed (1,1) double            {mustBeReal}               = NaN
        args.w      (1,1) struct                                       = struct()
        args.h      (1,1) struct                                       = struct()
    end

    x      = args.x;
    t_res  = args.t_res;
    f      = args.f;
    slot   = args.slot;
    x_feed = args.x_feed;

    if (isempty(x) || isempty(t_res) || isempty(f) || isempty(slot) || isnan(x_feed)) ...
            && (isempty(fieldnames(args.w)) || isempty(fieldnames(args.h)))
        eid = 'nonValidArguments:eval_vx';
        msg = 'The x, t_res, f, slot, and d_feed or w and h arguments must be specified.';
        
        error(eid, msg);
    end

    Nt = length(i);

    if ~isempty(x)
        dx = x(2) - x(1);
        Nx = length(x);
    else
        Nx = size(args.h, 1);
    end

    if ~isempty(f)
        Nf = length(f);
    end

    if isempty(fieldnames(args.w))
        W = (1 ./ eval_Zm(f, slot, 0, 'dx', dx)) .^ 2;
        w = 2 * real(eval_IFT(t_res, f, W));
    else
        W = args.w.W;
        w = args.w.w;
    end
    
    if isempty(fieldnames(args.h))
        Zx = NaN(Nx, Nf);
        for idx = 1 : Nx
            Zx(idx, :) = eval_Zm(f, slot, x(idx) - x_feed, 'dx', dx);
        end
    
        H = Zx .* W;
        
        h = NaN(Nx, Nt);
        parfor idx = 1 : Nx
            h(idx, :) = 2 * real(eval_IFT(t_res, f, H(idx, :)));
        end
    else
        h = args.h.h;
    end

    vx = NaN(Nx, Nt);
    
    % First time index
    vx(:, 1) = sum(i(:, 1) .* fliplr(h(:, 1)), 2) ./ w(:, 1);
    
    for m = 2 : Nt
        conv_part_1 = sum(i(:, 1 : m) .* fliplr(h(:, 1 : m)), 2);
        conv_part_2 = sum(vx(:, 1 : m - 1) .* fliplr(w(:, 2 : m)), 2);

        vx(:, m)    = (conv_part_1 - conv_part_2) ./ w(:, 1);
    end

    if isempty(fieldnames(args.w)) || isempty(fieldnames(args.h))
        varargout = {vx, struct('H', H, 'h', h), struct('W', W, 'w', w)};
    else
        varargout = {vx};
    end
end

