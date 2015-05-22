====== Projet 4: Prédiction de la structure secndaire ======

===== Introduction =====


La fonction d'une protéine et les intéractions entre plusieurs protéines dépendent directement de leur structure en trois dimensions qui elle-même dépend des différentes structures de segments locaux dans les protéines, appelées structure secondaire. Celle-ci est déterminée par les liens hydrogènes formés entre deux acides aminés (ou plus) de la protéine qui donnent au segment local une certaine structure. Les structures les plus courantes sont les α-helices, les ß-feuillet, les ß-tours et les bobines. 

Un des objectifs les plus difficiles et les plus importants de la bioinformatique est la prédiction de ces structures secondaires à partir d'une séquence de protéines afin de déterminer les fonctions et les liaisons possibles d'une protéine nouvellement séquencée.

L'algorithme GOR proposé par Jean Garnier et al. en 1978 analyse la probabilité qu'un acide aminé dans une séquence se trouve dans une certaine structure en fonction des acides aminés dans son entourage proche. Cet algorithme a été plusieurs fois amélioré depuis sa création (la dernière version en date est GOR V). J'ai implémenté ici la version III specifiée dans l'énoncé.

Le programme fonctionne en trois phases. Premièrement, il s'agit de parser un ensemble de fichier DSSP afin de récupérer un ensemble de séquences associées à leur structure. Deuxièment, le programme extrait des fréquences de cet ensemble de séquences. Pour terminer, le programme tente de prédire la structure d'un ensemble de séquences à partir de ces fréquences.

===== Les phases du programme =====

==== Parsing des fichiers DSSP ====

Dans un premier lieu on va devoir récupérer toutes les informations nécéssaires au bon fonctionnement de notre programme dans les fichiers fournis dans le cours.
On va utiliser le fichier "CATH_info.txt" afin de savoir quels fichiers dssp doivent être traités.

De plus le dernier caractère du nom du fichier dssp représente l'acide aminé dont il faut extraire la séquence (ex: 1B3A.dssp désigne le fichier 1B3.dssp dont il faut extraire les acides aminés de la chaîne A) car oui en effet il peut y avoir plusieurs chaînes de la même séquence, plusieurs structures d'une même protéine.

Ensuite après avoir déterminé quels fichiers dssp utiliser il faut parser l'information de ces dossiers qui nous intéresse:
  * le nom de la proteine
  * l'organisme dont il provient
  * la chaîne à laquelle appartient l'acide aminé
  * l'acide aminé
  * la structure secondaire de l'acide aminé

Pour les deux premiers, il suffit de lire la 4ème et 5ème, de split les éléments et de garder la partie intéressante.
Pour les 3 derniers si on essaye de split, on aura des problèmes car un des résultats attendus est un espace, or si on split l'espace ne sera pas pris en compte, il faut donc trouver une altérnative qui consiste en parser directement les cases de la ligne qui nous intéressent càd la 11ème case pour la chaîne, la 13ème case pour l'acide aminé et finalement la 16ème case pour la structure secondaire.

On va ensuite sauvegarder ces données dans un fichier avec la structure suivante:
> fichierParsé | nomProteine | organisme
sequence
structure secondaire

==== Les fréquences ====

=== Calcul des fréquences ===

Le programme parcoure le fichier créé lors de la phase précédente et à partir de chaque association séquence/structure, enregistre deux types de fréquences : la fréquence individuelle et la fréquence locale. Celles-ci seront cruciales lors de la troisième phase.

== Fréquence individuelle ==

La fréquence individuelle $$sfreq_{RS}$$ est le nombre de fois où le résidu R est impliqué dans une structure S. Elle permet de déterminer le score individuel, c'est-à-dire la probabilité que le résidu R se trouve dans une structure S.

== Fréquence locale ==

La fréquence locale $$pfreq_{RSmR_{m}}$$ est le nombre de fois où le résidu R est impliqué dans une structure S et situé à une distance m du résidu R_{m}. Cela permet de déterminer le score local, c'est-à-dire la probabilité que le résidu R se trouve dans une structure S en fonction des résidus de son entourage. Dans l'algorithme GORIII, les valeurs de m sont prises dans une fenêtre [-8,8] qui a été calculée à partir de la longueur moyenne des structures secondaires des protéines connues.

==== Prédiction de la structure ====

La troisième étape est la plus importante car elle effectue la prédiction. Le principe de base est de parcourir la séquence dont on veut prédire la structure et pour chaque résidu de calculer l'information d'appartenance à chaque structure. La structure dont l'information est la plus haute pour ce résidu sera la structure prédite.
Mon algorithme additionne deux informations pour obtenir cette information d'appartenance : l'information directionnelle et l'information de pair. Elles ont pour point commun l'utilisation de la "window". Pour chaque résidu, l'information dépendra donc des résidus qui l'entourent.
Ces informations sont calculées à partir du score individuel d'un résidu et du score local d'un résidu, eux-mêmes calculés à partir des fréquences obtenues lors de la phase précédente. Pour plus de détails, je vous renvoie à l'article de Jean Garnier et al. 1996 (http://www.ulb.ac.be/di/map/tlenaert/Home_Tom_Lenaerts/INFO-F-208_files/1996%20Garnier.pdf)

Dans les formules suivantes, $$I(S_{j};R_{j})$$ est le score individuel d'un résidu R à l'emplacement j pour une structure S et $$I(S_{j};R_{j+m}|R_{j})$$ est le score local d'un résidu R à l'emplacement j pour une structure S en tant compte d'une window [-m;+m]

== Information directionnelle ==

L'information directionnelle comprend le score individuel du résidu auquel est ajouté le score individuel des résidus dans la "window" pondéré par la distance à laquelle ils se trouvent du résidu.

$$I(S_{j};R_{1}..R_{m*2}) = I(S_{j};R_{j}) + \sum{m,m<>0}I(S_{m};R_{m})/abs(m)$$

== Information de pair ==
L'information de pair comprend le score individuel du résidu auquel est ajouté le score local des résidus dans la "window".

$$I(S_{j};R_{1}..R_{m*2}) = I(S_{j};R_{j}) + \sum{m,m<>0}I(S_{j};R_{j+m}|R_{j})$$


====Création du dictionnaire ====

Pour l'implémentation de GOR III, on a besoin de garder un certain nombre de couple de valeur pour le calcul des logarithmes. 
En effet, lorsqu'on parcourt une séquence par bloc afin d'incrémenter les valeurs de triples f(Structure, acide j, acide j + m), chaque valeur doit rester disponible lorsqu'on passe à l'indice suivant, à la séquence suivante.

Pour ce faire, on déclare un dictionnaire auquel va venir s'ajouter tous les couples f(Structure, acide j) lorsque les indices j et m sont égaux (lorsqu'on arrive à l'acide aminé au centre de la fenêtre d'évaluation) et tous les triples (Structure, acide j, acide j + m) lorsque j ≠ m.
Ces couples/triples de valeurs feront office de clefs dont les valeurs représenteront les occurences de ces combinaisons. 

Une fois toutes les séquences traitées, les combinaisons présentes dans le dictionnaire vont permettre de prédire la structure secondaire d'un acide aminé en particulier. Ce dictionnaire est envoyé à l'algorithme de calcul.

==== Algorithme GORIII ====

Pour chaque acide aminé d'une séquence, on peut calculer "l'information mutuelle" I(delta Sj ; Rj ; Rj+m) pour les 4 structures (Hélice alpha, beta brin, beta feuille, coude) grâce au dictionnaire. Cela mènera à quatre résultats différents et le plus grand des quatre est celui dont la conformation est la plus probable.

Donc en résumé, pour chaque acide de la séquence, on calcule son information mutuelle pour chaque structure possible. Le maximum des quatre représente la structure ayant le plus de potentiel pour être sa réelle structure secondaire. Une fois les prédictions faites sur tous les acides aminés, on peut renvoyer la séquence, sa structure secondaire réelle et la prédiction faite pour effectuer une analyse de la qualité des prédictions. 

==== Qualité des prédictions ====

A la fin du programme, il effectuer une mesure de qualité, c'est-à-dire que nous devons vérifier si les prédictions établies correspondent à la structure secondaire de la protéine en question. Deux calculs peuvent être effectués pour cela : la mesure Q3 et le coefficient de corrélation de Mathews.

  * Mesure Q3 : Un simple calcul entre le rapport du nombre de résidus correctements prédits et le nombre de résidus au total. Le Q3 n'est pas spécialement un bon indicateur car deux prédictions d'une même séquence peuvent avoir le même Q3 pour des raisons différentes : une prédiction (P1) peut s'être trompée sur les structures secondaires au bord d'un groupement spécifique mais avoir un bon nombre de structures secondaires. La prédiction (P2) aurait par exemple le même nombre de structures secondaires correctes que P1 mais les fautes se situraient en blocs (tous les betha-brins pouvant être identifiés comme des hélices-alpha par exemple). C'est pour cela qu'une seconde méthode pour mesurer la qualité des prédictions est nécessaire

  * MCC : Cette mesure de qualité va quant à elle s'appliquer sur la notion de vrai/faux négatif/positif concernant une structure bien précise. Le vrai positif indique que la structure secondaire de l'acide aminé en question a été correctement prédite. Le vrai négatif indique que la structure secondaire des autres acides aminés a été correctement prédite. Un faux positif indique que on a prédit la structure secondaire alors qu'elle n'apparait pas dans la structure réelle et un faux négatif indique que l'on a prédit une autre structure secondaire que celle prévue dans la structure réelle. Une fois chacune de ces données rassemblées, on peut calculer le MCC selon la formule suivante: $MCC = (TP * TN - FP * FN)/((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))^1/2$ avec TP = TruePositive, TN = TrueNegative, FP = FalsePositive and FN = FalseNegative.
Après avoir calculé chaque mesure de qualité pour chaque séquence, on peut calculer la qualité du modèle en temps que tel (en effectuant la qualité moyenne sur toutes).

==== Tests de l'algorithme ====

Alors on va tout d'abord parser et traiter "CATH_info.txt" qui seront nos données d'entraînement, ensuite on va faire de même avec "CATH_info_test.txt" qui seront nos données de test. Dans "CATH_info.txt"et "CATH_info_test.txt" on va récupérer la liste des fichier ".dssp" à parser et les chaînes qui devront être parsées dans les fichiers.
On utilise les données d'entraînement pour créer le dictionnaire qui va contenir les fréquences nécessaires pour le calcul des prédictions.


====Bibliographie (simplifiée)====
  * http://plage-desinvolte.pagesperso-orange.fr/d_agora/d_bioinfo/N-bioinfo.pdf
  * http://books.google.be/books?id=jhbkPHgRnj8C&pg=PA123&lpg=PA123&dq=gor+3+bioinformatique&source=bl&ots=uAyS_2-1hd&sig=9WNpw_zDLBo4RjLLx-Mna1JqE8U&hl=fr&sa=X&ei=y0tvU5X_A6_y7AbkiYGAAg&ved=0CFwQ6AEwBg#v=onepage&q=gor%203%20bioinformatique&f=false


