function Z = eval_Z(f, slot, d_feed)
%EVAL_Z Summary of this function goes here
%   Detailed explanation goes here
    arguments
        f      (1,:) double {mustBePositive}
        slot         SlotInDielectrics
        d_feed (1,1) double {mustBeReal}
    end
    
    Zin = slot.compute_zin(f);
    Zm  = eval_Zm(f, slot, d_feed);

    Nf = length(f);
    Z  = NaN(2, 2, Nf);
    
    Z(1, 1, :) = Zin;
    Z(1, 2, :) = Zm;
    Z(2, 1, :) = Zm;
    Z(2, 2, :) = Zin;
end
