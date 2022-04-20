% normalisation function
function matricenorm = normalisation(matrice)
maxi = max(max(matrice));
mini = min(min(matrice));
matricenorm = (matrice-mini)/(maxi-mini);