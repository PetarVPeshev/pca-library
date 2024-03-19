function Z = eval_Z(f, slot, d_feed, args)
%EVAL_Z Summary of this function goes here
%   Detailed explanation goes here
    arguments
        f                 (1,:) double             {mustBePositive}
        slot                    SlotInDielectrics
        d_feed            (1,:) double             {mustBeReal}
        args.MaxNumerical (1,1) double             {mustBeNonnegative} = 0
    end
    
    Nf = length(f);
    Nd = length(d_feed);

    Zin = slot.compute_zin(f);
    Z   = NaN(2, 2, Nf, Nd);

    for idx = 1 : Nd
        if d_feed(idx) <= args.MaxNumerical || args.MaxNumerical == 0
            Zm = eval_Zm(f, slot, d_feed(idx));
        else
            Zm = eval_Zmp(f, slot, d_feed(idx));
        end
        
        Z(1, 1, :, idx) = Zin;
        Z(1, 2, :, idx) = Zm;
        Z(2, 1, :, idx) = Zm;
        Z(2, 2, :, idx) = Zin;
    end
end

function Z = create_Z_matrix(f, slot, d_feed, args)
%EVAL_Z Summary of this function goes here
%   Detailed explanation goes here
    arguments
        f                 (1,:) double             {mustBePositive}
        slot                    SlotInDielectrics
        d_feed            (1,1) double             {mustBeReal}
        args.NumFeeds     (1,1) double             {mustBeReal, mustBePositive} = 2
        args.MaxNumerical (1,1) double             {mustBeNonnegative}          = 0
    end
    
    Nf = length(f);

    z_self = slot.compute_zin(f);
    Z   = NaN(2, 2, Nf, Nd);

    for idx = 1 : Nd
        if d_feed(idx) <= args.MaxNumerical || args.MaxNumerical == 0
            Zm = eval_Zm(f, slot, d_feed(idx));
        else
            Zm = eval_Zmp(f, slot, d_feed(idx));
        end
        
        Z(1, 1, :, idx) = Zin;
        Z(1, 2, :, idx) = Zm;
        Z(2, 1, :, idx) = Zm;
        Z(2, 2, :, idx) = Zin;
    end
end
