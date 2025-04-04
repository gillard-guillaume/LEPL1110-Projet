# Finite Element Mesh Generation and Finite Element Solver in C

## Auteurs

Projet réalisé dans le cadre du cours **LEPL1110**.

- **Gillard Guillaume**
- **Grégoire Igor**

Encadré par :
- **Tuteur** : Guilmain Vithit
- **Assistants** : Quiriny, Antoine 
- **Professeur responsable** : Vincent Legat


## Introduction

Ce projet a pour objectif d'implémenter une simulation la déformation d'une aile générée par transformation de Joukovsky face à des forces extérieurs en utilisant la méthode des éléments finis, dans le cadre du cours **LEPL1110**.  
L'aile est tout d'abord créée à partir d'une transformation conforme de Joukowsky, puis un maillage est généré à l'aide du logiciel **Gmsh**.  
Enfin, les forces sont appliquées et la déformation de l'aile est simulée en utilisant la méthode des éléments finis développée en **C**.


> ⚠️ **Attention** : Le SDK de Gmsh intégré fonctionne uniquement sous Linux.  
> Pour d'autres systèmes d'exploitation (Windows, macOS), veuillez télécharger et installer Gmsh depuis le site officiel : [https://gmsh.info/](https://gmsh.info/)

---
## Utilisation
  

Modifier les paramètres dans **src/main**
- Profil de Joukovsky
  - `R`: Rayon du cercle
  - `mu_x` : Centre du cercle sur l'axe des réels
  - `mu_y` : Centre du cercle sur l'axe des imaginaires
- Paramètres du matériaux (par défaut Aluminium 7075-T6)
  - `E` : Module de Young
  - `nu` : Module de Poisson
  - `rho` : Masse volumique
- Paramètre d'environnement
  - `g` : Gravitation terrestre
  - `g_multiplier` : Multiplicateur de gravité
- Paramètres du problème 
  - `renumType` : WIP
  - `solverType` : Type de solveur parmis **Cholvesky** (FEM_CHOV), **Gauss** (FEM_GAUSS) et **gradient conjugué** FEM_CG
  - `caseType` : Déformations planes ou tensions planes
  - `elementType` : Forme des éléments (**FEM_TRIANGLE**)

Pour lancer le programme :
```bash
# Créer un dossier de compilation
mkdir build
cd build

# Générer les fichiers de Makefile
cmake ..

# Compiler et exécuter
make run
```


---
## Arborescence


### 📁 fem/
- **fem.c / fem.h** :  
  Définition des structures principales (`femGeo`, `femDomain`, `femProblem`) et des fonctions d'assemblage du système éléments finis.

### 📁 forces/
- **forces.c / forces.h** :  
Implémentation du profil de force utilisé pour appliquer les conditions de Neumann normales sur les frontières du domaine.


### 📁 gl/
- **glfem.c / glfem.h** :  
  Fonctions de visualisation graphique utilisant OpenGL et GLFW : affichage du maillage, des déplacements et des déformations.

### 📁 mesh/
- **generateMesh.c / generateMesh.h** :  
  Fichier principal de la génération de maillage.  
  Il regroupe et appelle toutes les fonctions nécessaires pour créer la géométrie et générer le maillage. Cela permet de générer l'ensemble du maillage en appelant une seule fonction.

- **joukowsky.c / joukowsky.h** :  
  Construction du profil aérodynamique basé sur la transformation conforme de Joukowsky.

- **mesh.c / mesh.h** :  
  Définition des structures de base du maillage (nœuds, éléments, frontières) et outils de gestion.  
  Utilisation du moteur géométrique OpenCASCADE via Gmsh pour la génération de la géométrie et du maillage.

- **fixmesh.py** :  
  Script Python fourni par les assistants pour corriger et nettoyer les maillages générés.


### 📁 solver/
- **renum.c** :  
  Rénumération des nœuds pour améliorer les performances de résolution (minimisation du fill-in).
- **solver.c** :  
  Résolution du système linéaire issu de l'assemblage par un solveur direct.

### 📁 utils/
- **messages.c** :  
  Gestion de l'affichage des messages console (erreurs, informations, warnings).

### Fichiers principaux
- **main.c / main.h** :  
  Programme principal : gestion de l'enchaînement des étapes génération → assemblage → résolution → affichage.

---

## Prérequis

- Gmsh SDK (intégré, Linux uniquement)
- OpenGL
- GLFW
- BLAS / LAPACK (bibliothèques de calcul numérique)

---

## Installation

```bash
# Cloner le projet
git clone <git@github.com:gillard-guillaume/LEPL1110-Projet.git>

# Créer un dossier de compilation
mkdir build
cd build

# Générer les fichiers de Makefile
cmake ..

# Compiler et exécuter
make run
```

## Contrôles clavier

Lors de l'exécution du programme, différentes touches permettent d'interagir avec l'affichage :

- **`V`** : Afficher le maillage avec les déplacements (mode par défaut).
- **`S`** : Afficher la matrice associée au problème.
- **`D`** : Afficher uniquement le domaine géométrique (sans déformation).
- **`N`** : Changer le domaine affiché (uniquement en mode domaine (`D`)).
