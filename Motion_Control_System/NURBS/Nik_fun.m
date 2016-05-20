function Nik = Nik_fun(u , ui , uik1 , Nik1 , uik , ui1 , Ni1k1)
    Nik = 0;
    if  (uik1 - ui) ~= 0 && (uik - ui1) ~= 0
        Nik = (u - ui) / (uik1 - ui) * Nik1 + (uik - u) / (uik - ui1) * Ni1k1;
    elseif (uik1 - ui) ~= 0
        Nik = (u - ui) / (uik1 - ui) * Nik1;
    elseif (uik - ui1) ~= 0
        Nik = (uik - u) / (uik - ui1) * Ni1k1;
    end
end