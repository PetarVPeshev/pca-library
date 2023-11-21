function [magn_symb, magn_order] = find_magnitude(num)
    %FIND_MAGNITUDE Summary of this method goes here
    symbs = [ "Q", "R", "Y", "Z", "E", "P", "T", "G", "M", "k", "", ...
              "m", "u", "n", "p", "f", "a", "z", "y", "r", "q" ];
    
    magn = 1e30;
    magn_idx = 1;

    while magn > 1e-30
        if floor(num / magn) >= 1 && ceil(num / magn) < 1000
            magn_symb = convertStringsToChars(symbs(magn_idx));
            magn_order = magn;
            break;
        end

        magn = magn / 1e3;
        magn_idx = magn_idx + 1;
    end
end
