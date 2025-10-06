R Notebook
================

## **Installation du package “DADA2”**

``` r
library(dada2)
```

## **Chargement du package de données à traiter**

``` r
path <- "~/TUTORIEL-ADM/MiSeq_SOP"
list.files(path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

## **Création des listes des fichiers Forward (fnFs) et Reverse (fnRs)**

``` r
fnFs = sort(# Les met dans l'ordre alphabétique
  list.files(path, pattern="_R1_001.fastq", # Cherche dans le dossier "path" tous les fichiers qui contiennent "_R1_001.fastq"
                    full.names = TRUE)) # Garde le chemin complet

fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
```

## **Nom de chaque échantillon extrait du nom du fichier**

``` r
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

## **Visualisation des profils de qualité de lecture**

### **Qualité de lecture Forward**

``` r
plotQualityProfile(fnFs[1:2])
```

![](DADA2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

C’est un diagnostic visuel pour voir la qualité des bases (nucléotides)
le long des lectures.

- Le fond gris représente la fréquence de chaque score de qualité à
  chaque position de la lecture

- La ligne verte correspond à la moyenne du score de qualité à chque
  position

- La ligne orange correspond au quartiles montrant la variabilité de la
  qualité

- La ligne rouge signifie la proportion des lectures qui atteignent
  cette position (Illumina = plate, les lectures ont la même longeur)

–\> Ici, bonne qualité de lecture forward

## **Qualité de lecture Reverse**

``` r
plotQualityProfile(fnRs[1:2])
```

![](DADA2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Contrairement à la qualité des lectures forward, la qualité des lectures
reverse sont beaucoup moindre. On observe une chute brutale de la
qualité (en vert) ainsi que sa variation (en orange). Ce phénomène est
normal avec Illumina et DADA2 à un algorithme assez robuste pour les
séquences de moindre qualité en intégrant des informations de qualité
dans son modèle d’erreur.

## **Filtrer et Rogner**

``` r
filtFs <- file.path(path, "filtered", # Crée un chemin vers un sous-dosier "filtered" à l'intérieur du dossier 'path"
                paste0(sample.names,"_F_filt.fastq.gz")) # Contruit le nom du fichier filtré pour chaque échantillon

filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#Attribution des noms des échantillons aux fichiers (fichier filtré lié à son échantillon)
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

Les fichiers ont maintenant des nom liés directement au nom de
l’échantillon

### **Filtrage et pré-traitement des séquences**

``` r
# Prend en entrée les fichiers bruts "FASTQ" de forward et reverse, leur applique des filtres de qualité et les sauvegarde dans de nouveaux fichiers

out <- filterAndTrim(fnFs, # Fichier brut Forward
                     filtFs, # Fichier Forward filtré de sortie
                     fnRs, # Fichier brut Reverse
                     filtRs, # Fichier Reverse filtré de sortie
                     truncLen=c(240,160), # Tronque les lectures Forward à 240 bases et Reverse à 160 bases
              maxN=0, # Elimine toute lecture contenant un "N" (base inconnu)
              maxEE=c(2,2), # Max d'erreurs attendues par lecture (2 pour Forward et 2 pour Reverse)
              truncQ=2, # Coupe une lecture dès qu'un score de qualité tombe en dessous de 2
              rm.phix=TRUE, # Supprime les lectures provenant du génome PhiX (contrôle pour Illumina)
              compress=TRUE, # Enregistre les fichiers en ".fastq.gz"
              multithread=FALSE
              )

head(out) # Affiche un tableau résumant le filtrage
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

Le tableau montre clairement le nombre lectures avant filtrage
(reads.in) et après filtrage (reads.out). Ici, on observe que la
majorité des lectures est conservé (seulement un petite partie filtrée)

## **Taux d’erreur de séquençage**

Les séquenceurs font parfois des erreurs, certaines bases lues ne sont
pas les vraies bases. DADA2 permet de distinguer les vraies séquences
biologiques des séquences erronées.

–\> DADA2 = modèle d’erreur PARAMETRIQUE

``` r
errF <- learnErrors(filtFs, multithread=TRUE) # Apprend un modèle statistique des erreurs à partir des lectures filtrées Forward
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE) # Apprend un modèle statistique des erreurs à partir des lectures filtrées Reverse
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

## **Visualisation des taux d’erreur**

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    ## Transformation introduced infinite values in continuous y-axis

![](DADA2_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

- Les points correspondent aux taux d’erreurs observés
- La ligne noire correspond au taux d’erreur estimé par l’algorithme
- La ligne rouge correspond au taux d’erreur attendu selon le Q score
  (taux d’erreur théorique selon le Q score).

On observe bien que les points suivent bien la ligne noire et que par
conséquent, les taux d’erreurs observés sont partiquement égaux aux taux
d’erreur estimé. De plus, on observe bien une diminution du taux
d’erreur plus le Q score augmente qui tend à suivre la ligne rouge.

## **Application de l’algo de DADA2**

``` r
dadaFs <- dada(filtFs, # Fichiers forward filtrés
               err=errF, # Modèle d'erreurs appris précédemment
               multithread=TRUE # Permet d'utiliser plusieurs coeurs en même temps pour aller plus vite
               )
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

Pour chaque lecture filtrée, l’algorithme va comparer la séquence
observée au modèle d’erreurs. Il corrige les lectures en tenant compte
des probabilités d’erreurs :

- Si une différence peut-être expliquée par une erreur de séquençage =
  Séquence corrigée

- Si une différence est trop grande pour être une simple erreur =
  Séquence considérée comme une vraie variante biologique (ASV)

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

Idem que précédemment mais avec les lectures filtrées de Reverse

## **Inspection**

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

Exemple : Pour le 1er échantillon, l’algorithme à identifier 128 vrais
variants dans 1979 séquences

## **Fusion des lectures appariées**

``` r
# Reconstruction des séquences complètes des amplicons
mergers <- mergePairs(dadaFs, # ASVs du brin forward
                      filtFs, # Fichiers FASTQ filtrés forward
                      dadaRs, # ASVs du brin reverse
                      filtRs, # Fichiers FASTQ filtrés reverse
                      verbose=TRUE # Affiche des détails pendant l'exécution
                      )
```

    ## 6540 paired-reads (in 107 unique pairings) successfully merged out of 6891 (in 197 pairings) input.

    ## 5028 paired-reads (in 101 unique pairings) successfully merged out of 5190 (in 157 pairings) input.

    ## 4986 paired-reads (in 81 unique pairings) successfully merged out of 5267 (in 166 pairings) input.

    ## 2595 paired-reads (in 52 unique pairings) successfully merged out of 2754 (in 108 pairings) input.

    ## 2553 paired-reads (in 60 unique pairings) successfully merged out of 2785 (in 119 pairings) input.

    ## 3646 paired-reads (in 55 unique pairings) successfully merged out of 4109 (in 157 pairings) input.

    ## 6079 paired-reads (in 81 unique pairings) successfully merged out of 6514 (in 198 pairings) input.

    ## 3968 paired-reads (in 91 unique pairings) successfully merged out of 4388 (in 187 pairings) input.

    ## 14233 paired-reads (in 143 unique pairings) successfully merged out of 15355 (in 352 pairings) input.

    ## 10528 paired-reads (in 120 unique pairings) successfully merged out of 11165 (in 278 pairings) input.

    ## 11154 paired-reads (in 137 unique pairings) successfully merged out of 11797 (in 298 pairings) input.

    ## 4349 paired-reads (in 85 unique pairings) successfully merged out of 4802 (in 179 pairings) input.

    ## 17431 paired-reads (in 153 unique pairings) successfully merged out of 17812 (in 272 pairings) input.

    ## 5850 paired-reads (in 81 unique pairings) successfully merged out of 6095 (in 159 pairings) input.

    ## 3716 paired-reads (in 86 unique pairings) successfully merged out of 3894 (in 147 pairings) input.

    ## 6865 paired-reads (in 99 unique pairings) successfully merged out of 7191 (in 187 pairings) input.

    ## 4426 paired-reads (in 67 unique pairings) successfully merged out of 4603 (in 127 pairings) input.

    ## 4576 paired-reads (in 101 unique pairings) successfully merged out of 4739 (in 174 pairings) input.

    ## 6092 paired-reads (in 109 unique pairings) successfully merged out of 6315 (in 173 pairings) input.

    ## 4269 paired-reads (in 20 unique pairings) successfully merged out of 4281 (in 28 pairings) input.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                       sequence
    ## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
    ## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
    ## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
    ## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
    ## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       579       1       1    148         0      0      1   TRUE
    ## 2       470       2       2    148         0      0      2   TRUE
    ## 3       449       3       4    148         0      0      1   TRUE
    ## 4       430       4       3    148         0      0      2   TRUE
    ## 5       345       5       6    148         0      0      1   TRUE
    ## 6       282       6       5    148         0      0      2   TRUE

Permet d’observer les vraies séquences les plus abondantes dans
l’échantillon. Les colonnes “nmismatch” (nombre de bases différentes
dans la zone de chevauchement) et “nindel” (insertion ou délétions
observées) permettent de voir si la fusion est cleen.

## **Construction d’une table de séquence**

``` r
seqtab <- makeSequenceTable(mergers) 
dim(seqtab)
```

    ## [1]  20 293

Permet de créer une table de séquence qui prend en entrée les séquences
fusionnées.

- Lignes : Echantillons

- Colonnes : ASVs

- Valeurs : Abondance des séquences

La commande “dim()” permet de retourner la matrice de longueur 2. Ici,
le 1er chiffre correspond au nombre d’échantillon (20) et le 2ème
chiffre au nombre d’ASV (293)

### **Inspection de la distribution des longueurs des séquences**

``` r
table( # Compte combien de séquences ont une certaine longueur
  nchar( # Calcule la longueur de chaque séquence
    getSequences(seqtab) # Extrait tous les ASVs présentes dans "seqtab"
        ))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  88 196   6   2

Ici, on oberve que :

- 1 séquence à une longueur de 251 nucléotides

- 88 de 252 nucléotides

- 196 de 253 nucléotides

- 6 de 254 nucléotides

- 2 de 255 nucléotides

## **Supprimer les chimères**

Une chimère est une séquence artificielle formée quand 2 fragments d’ADN
s’hybrident partiellement et sont amplifiés comme une “fausse” séquence.

``` r
seqtab.nochim <- removeBimeraDenovo( #Enlève les séquences chimériques
  seqtab, # Table de séquence fusionnée
  method="consensus", # Compare les séquences échantillon par échantillon et garde celles qui sont considérées "non-chimérique" dans la majorité des cas
  multithread=TRUE, # Accélère le calcul
  verbose=TRUE # Affiche un résumé
  ) 
```

    ## Identified 61 bimeras out of 293 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  20 232

Ici, 61 séquences sur les 293 ont été identifié comme des chimères

## **Proportion de lecture “survivante” au filtrage de chimères**

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.9640374

On observe que 96 % des lectures ont été conservé et que 4% sont des
chimères

## **Suivre les lectures à travers la pipeline**

``` r
getN <- function(x) sum(getUniques(x)) #Nombre de lectures uniques conservées pour 1 échantillon

track <- cbind(out, # Lectures initiales + filtrées
               sapply(dadaFs, getN), # Lectures après débruitage forward
               sapply(dadaRs, getN), # Lectures après débruitage reverse
               sapply(mergers, getN), # Lectures après fusion des séquences
               rowSums(seqtab.nochim) # Lectures finales non chimériques
               )

# Nom pour les colonnes et lignes du tableau
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# Visualisation
head(track)
```

    ##        input filtered denoisedF denoisedR merged nonchim
    ## F3D0    7793     7113      6976      6979   6540    6528
    ## F3D1    5869     5299      5227      5239   5028    5017
    ## F3D141  5958     5463      5331      5357   4986    4863
    ## F3D142  3183     2914      2799      2830   2595    2521
    ## F3D143  3178     2941      2822      2868   2553    2519
    ## F3D144  4827     4312      4151      4228   3646    3507

Permet une visualisation et traçabilité complète du pipeline DADA2. Il
montre combien de lectures ont garde/perds à chaque étape et par
échantillon. Normalement, on observe une petite perte à chaque étape ce
qui est normal (ici, c’est le cas)

## **Attribution d’une taxonomie**

Comparaison des ASVs à une base de données de références pour attribuer
une taxonomie

``` r
taxa <- assignTaxonomy(seqtab.nochim,  "~/TUTORIEL-ADM/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1", 
                       multithread=TRUE 
                       )
```

## **Examination des attributions**

``` r
taxa.print <- taxa # Copie l'objet "taxa" pour l'affichage
rownames(taxa.print) <- NULL # Supprime les noms de lignes (séquences d'ADN brutes)
head(taxa.print) # Affiche les 6 premières lignes
```

    ##      Kingdom    Phylum         Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus        
    ## [1,] NA           
    ## [2,] NA           
    ## [3,] NA           
    ## [4,] NA           
    ## [5,] "Bacteroides"
    ## [6,] NA

Les *Bacteroidota* sont très bien représentés parmi les taxons les plus
abondants.

## **Evaluation de la précision (échantillon de contrôle “Mock”)**

``` r
unqs.mock <- seqtab.nochim["Mock",] # Sélectionne uniquement la ligne correspondant à l'échantillon "Mock"

unqs.mock <- sort(unqs.mock[unqs.mock>0], # Ne garde que les ASVs présent
                  decreasing=TRUE # Trie les ASVs par abondance décroissante
                  ) 

cat("DADA2 inferred", length(unqs.mock), # Nombre d'ASVs différentes détectées dans l'échantillon "Mock"
    "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

Permet de vérifier la précision et la sensibilité de DADA2 sur le
contrôle “Mock” connu. Ici, DADA2 a identifié 20 séquences d’échantillon
présentes dans la communauté Mock.

## **Vérification des séquences détectées dans Mock**

``` r
# Extrait toutes les séquences d'ADN présentes dans le fichier "FASTA" de référence
mock.ref <- getSequences( file.path(path, "HMP_MOCK.v35.fasta"))

# Comparer les ASVs détectées aux séquences de référence
match.ref <- sum( # Combien de séquences détectées sont exactement présentes dans le fichier de référence
  sapply(
  names(unqs.mock), # Séquences des ASVs détectées dans "Mock"
  function(x) any( # Renvoie TRUE si au moins 1 correspondance est trouvée
    grepl(x, mock.ref) # Vérifie si cette séquence existe exactement dans les séquences de référence
                              )))

cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

20 séquences détectées par DADA2 et 20 séquences correspondent
exactement aux séquences connues de “Mock”. Cela signifie que le taux
d’erreur est de 0 % et qu’il n’y a pas d’erreur ni contaminant dans le
contrôle et dans les autres échantilllons.
