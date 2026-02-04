# Projet-C-GLV

Implémentation en C d'opérations sur les courbes elliptiques et de l'accélération GLV (Gallant–Lambert–Vanstone) pour la multiplication scalaire.

## Contenu du dépôt

- **Opérations de base sur courbes elliptiques** : addition en coordonnées affines et projectives, doublement, etc.
- **Multiplication scalaire** : méthodes classique (square-and-multiply) et accélérée via GLV.
- **Décomposition GLV** : calcul des scalaires courts et endomorphismes associés.
- **Outils annexes** : pré-calculs, solveur quadratique, vecteurs courts, génération de constantes (scripts Sage).

## Structure (fichiers principaux)

- `main.c` : point d’entrée et démonstrations.
- `EC_*.c` / `EC_*.h` : structures et opérations sur courbes elliptiques.
- `glv_*.c` / `glv_*.h` : décomposition et accélération GLV.
- `double_scalar_multiplication.c` : multiplication double.
- `precompute_table.c` : pré-calculs pour accélération.
- `gen_const_example*.sage` : génération de constantes.

## Compilation

Ce projet utilise un **Makefile**.

```bash
make
```

## Exécution

Après compilation, lancez l’exécutable généré 'glv' :

```bash
./glv
```

## Notes

- Le dossier nécessite la bibliothèques GMP (mpz_t). Vérifiez le Makefile si nécessaire.
- Les fichiers `.sage` servent à générer des constantes utilisées dans le code C.