# REQUIREMENTS.md — LabStats GUI

## Outil d'analyse statistique pour le labo (Windows, usage étudiant)

---

## 1. ARCHITECTURE GÉNÉRALE

- **Frontend (GUI) :** Python (Tkinter)
- **Backend (Statistiques) :** R via `subprocess` appelant `run_analysis.R`
- **Distribution :** PyInstaller (.exe Windows) + `run_analysis.R` dans le même dossier
- **R doit être installé** sur la machine. Le GUI doit détecter le chemin de `Rscript` automatiquement (PATH Windows) ou permettre à l'utilisateur de le spécifier une seule fois (sauvegardé dans un fichier de config local).

### Packages R requis (à documenter dans un README séparé pour l'installation)

- `readxl` (lecture Excel)
- `glmmTMB` (GLMM : NB, Gamma, Beta, ZOIB)
- `DHARMa` (validation des résidus)
- `emmeans` (moyennes marginales estimées + contrastes)
- `performance` (R² de Nakagawa)
- `ggplot2` + `ggeffects` (visualisation publication-ready)
- `knitr` / `rmarkdown` ou `officer` (génération PDF/rapport)
- `openxlsx` ou `readxl` (lecture .xlsx)
- `dplyr` (manipulation données)

---

## 2. INTERFACE GRAPHIQUE (GUI Python)

### 2.1 Chargement des données

- **Bouton "Parcourir"** : sélecteur de fichier `.xlsx` (Excel).
- **Dropdown "Feuille (Sheet)"** : liste les feuilles disponibles dans le fichier Excel sélectionné. L'utilisateur choisit la bonne feuille.
- Une fois la feuille sélectionnée, le GUI lit les noms de colonnes et peuple tous les dropdowns suivants.

### 2.2 Structure expérimentale

- **Dropdown "Animal ID"** : colonne contenant l'identifiant unique de chaque sujet (random intercept obligatoire).

### 2.3 Variables Indépendantes (VI) — Effets Fixes

- **3 Dropdowns "VI 1 / VI 2 / VI 3"** : sélection des colonnes de variables indépendantes.
- Chaque VI a un **label de type** (automatique ou manuel) : `Between-subject` ou `Within-subject`.
  - **Between-subject** (ex: Sexe, Phénotype, Génotype) : ne varient PAS à l'intérieur d'un animal → **Effet fixe seulement, pas de Random Slope**.
  - **Within-subject** (ex: Drug Treatment, Day, Session) : varient à l'intérieur d'un animal → **Candidats au Random Slope**.
- **Note importante pour le script R :** Seules les VI de type `Within-subject` seront tentées comme Random Slopes dans la structure aléatoire.

### 2.4 Variables Dépendantes (VD)

- **Jusqu'à 10 lignes de VD**, chacune composée de :
  - **Dropdown "VD"** : sélection de la colonne.
  - **Bouton radio / Dropdown "Type"** : `Counts` | `RT` | `Proportion`.
- Le type détermine la famille de distribution utilisée dans R (voir Section 4).

### 2.5 Actions

- **Bouton "Analyser"** : lance l'analyse (appelle `Rscript run_analysis.R` avec les paramètres).
- **Bouton "Rapport PDF"** : ouvre le PDF généré automatiquement dans le dossier `outputs/`.
- Une **barre de progression** ou un **message de statut** doit indiquer que R est en train de tourner (ex: "Analyse en cours... veuillez patienter").
- En cas d'erreur R, afficher un message clair dans le GUI (pas juste un crash silencieux).

---

## 3. COMMUNICATION PYTHON → R

### 3.1 Mécanisme

Python appelle R via `subprocess.run()` :

```
Rscript run_analysis.R --input "chemin/fichier.xlsx" --sheet "Feuille1" --id "animal_id" --vi1 "drug" --vi1_type "within" --vi2 "day" --vi2_type "within" --vi3 "sex" --vi3_type "between" --vd1 "lever_presses" --vd1_type "counts" --vd2 "latency" --vd2_type "rt" --output_dir "outputs/"
```

### 3.2 Fichiers de sortie générés par R

- `outputs/report_YYYYMMDD_HHMMSS.pdf` : rapport complet.
- `outputs/results_summary.json` : résumé des résultats (pour affichage éventuel dans le GUI).
- `outputs/figures/` : dossier contenant les figures `.png` intégrées dans le PDF.

---

## 4. LOGIQUE STATISTIQUE (Script R — run_analysis.R)

### 4.1 Arbre décisionnel par type de VD

#### TYPE "Counts" (Comptes entiers : 0, 1, 2, 3...)

**Famille cible : GLMM Négative Binomiale (NB), lien Log**

1. Vérifier présence de zéros excessifs → test Zero-Inflation (DHARMa `testZeroInflation()`).
   - Si Zero-Inflation détectée → passer à **Zero-Inflated NB (ZINB)** via `glmmTMB(..., ziformula = ~1)`.
2. Ajuster le modèle NB.
3. Vérifier surdispersion → si NB ne suffit pas, noter dans le rapport.
4. Validation : DHARMa QQ-plot + KS test.
   - Si KS échoue → tenter ZINB, sinon signaler dans le rapport.

#### TYPE "RT" (Temps de réaction / Latences : valeurs continues positives)

**Famille cible : GLMM Gamma, lien Log**

1. Vérifier absence de zéros (Gamma ne tolère pas les zéros).
   - Si zéros présents → avertissement dans le rapport + suggestion de nettoyage.
2. Ajuster le modèle Gamma.
3. Validation : DHARMa QQ-plot + KS test.
   - Si KS échoue → tenter LMM sur log(VD) comme alternative, noter dans le rapport.

#### TYPE "Proportion" (Pourcentages continus : 0.0 à 1.0)

**Famille cible : GLMM Beta, lien Logit**

1. Convertir automatiquement si les données sont en % (0-100) → diviser par 100.
2. Vérifier présence de 0 exacts ou 1 exacts.
   - Si présents → passer à **Zero-One Inflated Beta (ZOIB)** via `glmmTMB(..., family = beta_family(), ziformula = ~1)`.
   - Sinon → **Beta standard**.
3. Validation : DHARMa QQ-plot + KS test.
   - Si KS échoue → signaler dans le rapport avec suggestion.

### 4.2 Structure aléatoire (pour tous les types)

**Stratégie en entonnoir (Étage 3 + 4 de l'arbre décisionnel) :**

1. **Modèle Maximal :** `(1 + VI_within | Animal_ID)` — Intercept + Slopes pour toutes les VI Within-subject.
2. **Si Singular Fit :** Retirer les corrélations → `(1 + VI_within || Animal_ID)`.
3. **Si encore Singular Fit :** Retirer les Random Slopes → `(1 | Animal_ID)`.
4. **Si Failure to Converge :** Changer d'optimiseur (BFGS), augmenter les itérations.
5. **Si échec total :** Passer au **Plan B : GLM + Cluster-Robust SE** (sandwich estimator). Signaler clairement dans le rapport que le modèle mixte n'a pas convergé et que cette correction a été appliquée.

**Note :** Les VI de type `Between-subject` (ex: Sexe, Phénotype) sont incluses comme effets fixes uniquement. Elles ne sont jamais tentées comme Random Slopes.

### 4.3 Tailles d'effets (Étage 5)

- **R² de Nakagawa** (`performance::r2()`) : rapporter R²_marginal et R²_conditionnel.
- **Moyennes Marginales Estimées (EMMeans)** via `emmeans` : pour chaque contraste significatif.
- **Taille d'effet par contraste :**
  - Counts / RT (NB / Gamma) → **Incidence Rate Ratio (IRR)** ou **Response Ratio**.
  - Proportion (Beta) → **Odds Ratio**.
- **Corrections pour comparaisons multiples :** FDR (Benjamini-Hochberg) par défaut sur les contrastes post-hoc.

### 4.4 Validation DHARMa (pour tous les types)

Pour chaque modèle ajusté :

- Générer les résidus simulés (`simulateResiduals()`).
- QQ-plot des résidus → sauvegarder en `.png`.
- Test KS (`testUniformity()`).
- Test de surdispersion (`testDispersion()`).
- Test de Zero-Inflation (`testZeroInflation()`).
- Si un test échoue → noter le problème ET la solution recommandée dans le rapport.

---

## 5. RAPPORT PDF

### 5.1 Structure du rapport (dans l'ordre)

#### Section 0 — Paramètres de l'analyse (Reproductibilité)

Reproduire exactement ce qui a été entré dans le GUI :

- Fichier Excel + Feuille utilisée.
- Colonne Animal ID.
- VI 1, 2, 3 + leur type (Between/Within).
- VD analysées + leur type (Counts/RT/Proportion).
- Date et heure de l'analyse.
- Version de R + versions des packages utilisés.

#### Section 1 — Méthodes (Texte auto-généré)

Paragraphe de méthodes prêt à copier-coller dans un article scientifique, du type, mais en anglais :

> *"Les données de [VD] ont été analysées à l'aide d'un modèle mixte généralisé (GLMM) avec une distribution [Négative Binomiale / Gamma / Beta] et un lien [Log / Logit]. La structure aléatoire incluait [description de la structure retenue]. La validité du modèle a été confirmée par l'inspection des résidus quantiles simulés via le package DHARMa (test de Kolmogorov-Smirnov, p = X.XX). Les tailles d'effets sont rapportées sous forme de [IRR / Odds Ratio] avec intervalles de confiance à 95%. Les comparaisons multiples ont été corrigées par la méthode FDR (Benjamini-Hochberg)."*

#### Section 2 — Résultats par VD

Pour chaque VD analysée :

- **Tableau des effets fixes** : Estimate, SE, z/t, p-value, IC 95%.
- **Mise en évidence visuelle** des effets significatifs (p < 0.05 en gras, p < 0.001 marqué ***).
- **R² de Nakagawa** : R²_marginal et R²_conditionnel.
- **Tableau EMMeans** : Moyennes marginales estimées par groupe/condition + IC 95%.
- **Tailles d'effets** : IRR ou Odds Ratio par contraste.
- **Figure 1 :** QQ-plot DHARMa (validation du modèle).
- **Figure 2 :** Graphique des EMMeans avec barres d'erreur (IC 95%) — publication-ready (ggplot2, thème minimal).

#### Section 3 — Avertissements et Diagnostics

Si des problèmes ont été détectés, les lister clairement :

- ⚠️ "Singular Fit détecté : la structure aléatoire a été simplifiée de [modèle maximal] à [modèle retenu]."
- ⚠️ "Zero-Inflation détectée : un modèle ZINB/ZOIB a été utilisé à la place du modèle standard."
- ⚠️ "Le test KS de DHARMa est significatif (p = X.XX) : le choix de distribution [X] pourrait ne pas être optimal. Considérer [solution recommandée]."
- ⚠️ "Le modèle mixte n'a pas convergé : un GLM avec Cluster-Robust SE a été utilisé (Plan B)."
- ⚠️ "Des valeurs exactement à 0 ou 1 ont été détectées dans [VD] : un modèle ZOIB a été appliqué."

#### Section 4 — Comment rapporter ces résultats

Petite section pédagogique avec une phrase-exemple de résultats pour chaque VD :

> *"Le traitement a significativement augmenté le nombre d'appuis sur le levier (IRR = 1.45, IC 95% [1.12, 1.87], p = 0.005, R²_m = 0.18)."*

---

## 6. CONTRAINTES TECHNIQUES

- **Langue du rapport :** Anglais (tous les textes auto-générés en anglais).
- **Encodage :** UTF-8 partout (important pour le français, il peut y avoir des caractères francais dans les données).
- **Gestion des erreurs :** Toute erreur R doit être capturée et affichée de façon lisible dans le GUI (pas de crash silencieux).
- **Timeout :** Si R tourne plus de 5 minutes, afficher un avertissement dans le GUI.
- **Noms de colonnes :** Tolérer les espaces et accents dans les noms de colonnes Excel (fréquent en contexte francophone).
- **Données manquantes (NA) :** Le script R doit les gérer silencieusement (GLMM les gère nativement). Rapporter le nombre de NA par VD dans la Section 0 du rapport.

---

## 7. CE QUI N'EST PAS DANS CE SCOPE

- Pas d'analyse de puissance (power analysis) — hors scope pour l'instant.
- Pas de détection automatique du type de VD — l'étudiant choisit manuellement (Counts / RT / Proportion).
- Pas de gestion multi-ordinateurs / réseau — usage sur un seul poste Windows de labo.
- Pas de base de données — chaque analyse est indépendante.
