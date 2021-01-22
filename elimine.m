
function [tilde_AA,tilde_LL] = elimine (AA,LL,Refneu, Nbpt, Coorneu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elimine :
% Procédure de pseudo-élimination des noeuds i correspondant à Refneu(i)=1
%
% SYNOPSIS tilde_AA,tilde_LL = elimine (AA,LL,Refneu)
%          
% INPUT * AA,LL,Refneu : AA et LL vont subir la pseudo-élimination selon les valeurs
%                        de Refneu
%
% OUTPUT - tilde_AA,tilde_LL: matrices pseudo-éliminées
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbpt = size(Refneu)(1);

vecteur_nul = zeros(Nbpt,1);


% Pour le premier exercice.
for i=1:Nbpt  %on parcourt les sommets
  if Refneu(i,1)==1  %si la référence du noeud est 1, on applique la pseudo-élimination
    AA(i,:)=vecteur_nul;  %on élimine la ligne i de AA
    AA(:,i)=vecteur_nul;  %on élimine la colonne i de AA
    AA(i,i)=1;  %on place le coefficient diagonal i de AA à 1
    LL(i)=0;  %on élimine la i-ème composante de LL
  endif
endfor

% Pour le deuxième exercice, on n'elimine rien.

%{

% Pour le dernier exercice.
% On choisit l'arête à droit de Omega_1 comme la partie du bord sous condition
% de Fourier. Alors on elimine les points qui ne sont pas sur cette arete.

for i=1:Nbpt  % On parcourt les sommets
  % Si la référence du noeud est 1 est le noeud n'est pas sur l'arete choisi,
  % on applique la pseudo-elimination.
  if Refneu(i,1) == 1 && Coorneu(i, 1) != 2
    AA(i,:)=vecteur_nul;  %on élimine la ligne i de AA
    AA(:,i)=vecteur_nul;  %on élimine la colonne i de AA
    AA(i,i)=1;  %on place le coefficient diagonal i de AA à 1
    LL(i)=0;  %on élimine la i-ème composante de LL
  endif
endfor

%}

%tilde_AA,tilde_LL = AA,LL;
tilde_AA = AA;
tilde_LL = LL;

endfunction
