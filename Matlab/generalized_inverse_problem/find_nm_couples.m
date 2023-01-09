function [liste_couples,liste_wave] = find_nm_couples(k1,k2,p,q)
%Find all couples (n,m) non-negative where
%   k1^2 <= pi^2((n/p)^2+(m/q)^2)<k2^2
%On cherche n entre p*k1 et p*k2 et m entre q*k1 et q*k2
nbr_points=0;liste_couples=[];liste_wave=[];
for n=0:round(k2*p)+1
    for m=0:round(k2*q)+1
        if (k1^2 <= pi^2*((n/p)^2+(m/q)^2) && pi^2*((n/p)^2+(m/q)^2) <= k2^2)
            if (nbr_points==0)
                liste_couples = [n,m];nbr_points=nbr_points+1;
                liste_wave = [pi^2*(n^2/p^2+m^2/q^2)];
            else
                liste_couples(end+1,1:2) = [n,m];
                liste_wave(end+1) = [pi^2*(n^2/p^2+m^2/q^2)];
            end
        end
    end
end
end

