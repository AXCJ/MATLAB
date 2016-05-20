function Nik_1 = Nik_md_fun(ui , uik1 , Nik1 , uik , ui1 , Ni1k1 , k)
    Nik_1 = 0 * Nik1;
    if  (uik1 - ui) ~= 0 && (uik - ui1) ~= 0
        Nik_1 = (k-1) * (Nik1 / (uik1 - ui) - Ni1k1 / (uik - ui1));
    elseif (uik1 - ui) ~= 0
        Nik_1 = (k-1) / (uik1 - ui) * Nik1;
    elseif (uik - ui1) ~= 0
        Nik_1 = (1-k) / (uik - ui1) * Ni1k1;
    end
end