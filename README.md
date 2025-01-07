# Projet conductivité thermique d'une maison

## Configurations différentes
| Couleur toiture | Alpha | uS (°C) | f |
|-|-|-|-|
| Blanc | alpha0 = 2e-1 | 1.2 * uE = 40 | 0 |
| Noir | 10 * alpha0 | 1.3 * uE | 0 |
| Fenêtre | 10e6 * alpha0 | 1.1 * uE | 0.25 |

## Paramètres
### Alpha
Conductivité thérmique sur le bord du maillage (ex : mur en béton).

### Kappa
Conductivité thérmique dans le maillage (ex : air ou eau).

### uS
Température exterieur.

### f
Source de charleur dans le maillage (ex : radiateur ou soleil qui frappe le sol).

## Maillage
Création d’un maillage sous forme de maison avec 5 bords : 2 murs, 2 toits et un sol. <br>
Le pas et la qualité du maillage sont très corrects, il n’y a pas de triangle dégénéré. <bR>
![image](/images/millage.png)
- h=0.3092164366780396
- Q=1.8225720390498603

## Résultat
### Toiture blanche
![image](/images/toit_blanc.png)
- Température max : 47.61861966878903
- Température min : 27.405391526679534
- Température moyenne : 38.70904176087704

### Toiture noir
![image](/images/toit_noir.png)
- Température max : 52.029163358606
- Température min : 25.348703676020413
- Température moyenne : 41.43358363521184

### Toiture fenêtre
![image](/images/toit_vitre.png)
- Température max : 44.00000004406034
- Température min : 25.000000317283277
- Température moyenne : 37.27965073801681
