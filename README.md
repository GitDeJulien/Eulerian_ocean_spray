# Eulerian_ocean_spray

## Pré-réglage (pour le développement):
> [!NOTE]
> Pour que le code fonctionne bien sous vscode il faut aller dans les paramètres en bas à gauche de la fenêtre vscode. 
> Cliquer sur __Settings__ (ou __Ctr+__ directement).
> Taper __Fortran__ dans la barre de recherche. 
> Aller dans __User->Lintings__.
> Ajouter __${workspaceFolder}/include__ dans la section __Mod Output__. 
> Une fois cela fais s'il reste des fichiers _.mod_ dans src/ ou sub/ vous pouvez les suprimer.

# Déscription:
### data/data.toml
 - Fichier _.toml_ regroupant les données utiles

### sub/mod_toml_parser.f90
 - Module permettant de lire un fichier _.toml_

### sub/mod_precision.f90
 - Définition d'une précision double _pr = 8_
 - Définition de constant (_pi_ et _g_)

### src/mod_data.f90
 - Module définissant une structure stockant toute les données

### src/mod_structured_mesh.f90
 - Module définissant un maillage structuré. Il peut être 1D ou 2D en fonction des paramètres rentrés
 - Ce module permet aussi de définir une strucure de type maille groupant les caractéristique d'une maille
 mais aussi la solution au centre de cette maille

### src/functions.f90
 - Module permettant d'initialiser toutes les coefficients de l'équation de Williams-Boltzmann

### src/sample.f90
 - Module définissant le tirage celon $N(r)$ avec accéptation-rejet

### src/init.f90
 - Module permettant d'initialiser la solution dans chaque cellule du maillage

### src/flux.f90
 - Module définissant un flux décentré 1D pour le moment 

### main.f90
 - Fichier main 

# Taches
 - [ ] Écrire la boucle en temps en itérant sur chaque cellule du maillage
 - [ ] Écire un fichier de sortie pour la solution à mettre dans output (format ?)
 - [ ] Vérifier que le passage du Lagrangien au Eulerien (x_p -> x) se fait avec des moyennes
 - [ ] Vérifier que l'initialisation du vecteur solution est la bonne. J'en doute...
 - [ ] Coder les fonctions manquante pour le 2D s'il y a le temps...

# Compilation et exécution
 - Pour compiler deux possibilités:
  - En mode debug: __make debug__ ou simplement __make__
  - En mode released: __make release__
 - Pour exécuter deux possibilités:
  - En mode debug: __make exe_debug__
  - En mode released: __make exe_release__
 - Pour clean:
  - __make clean__

> [!NOTE]
> Le repertoire __build__ acueille tous les fichier _.o_ ainsi que l'executable et le repertoire 
> __include__ acueille tous les fichiers _.mod_. 

