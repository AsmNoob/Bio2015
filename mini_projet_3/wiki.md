======Projet 3: La construction de profils======

===== Introduction aux PSSM =====


Projet à télécharger: {{:f20815:start:333083:mini_projet_3.zip|}}

Une PSSM (Position Specific Scoring Matrix) permet de représenter un motif dans un groupe de séquences (de protéines dans notre cas). Elles sont très utiles pour représenter une famille ou un domaine de protéines dans son ensemble.

Ce troisième projet consiste en l'implémentation d'un algorithme de création d'une PSSM et dans la vérification de sa pertinence biologique. Le principe est semblable à la construction d'une matrice de substitution mais en calculant le score pour chaque position de chaque acide aminé dans le groupe de séquence donné. On travaille avec un alignement multiple de séquences provenant du même domaine ou de la même famille dans le but de créer une PSSM spécifique à ce domaine ou cette famille.

Pour commencer, je vais faire une introduction des deux familles utilisés, s'en suivra une présentation de mon algorithme et de mon implémentation.

===== Présentation des familles =====

==== SH2 Family ====

=== Introduction au domaine SH2 ===

Le premier domaine que j'ai analysé est le domaine SH2. Ce domaine est contenu principalement dans l'oncoprotéine Src dont il tire son nom Src Homologuous 2. Son rôle principal est la communication cellulaire mais on le trouve également dans des protéines impliquées dans la transcription ou dans des enzymes commet les phosphatases dans lesquelles son but est la détection tyrosines.

Le domaine SH2 permet à une protéine de s'attacher à d'autres protéines par un résidu tyrosine, c'est pourquoi on le trouve souvent dans des protéines impliquées dans la transduction d'un signal par les récepteurs tyrosines kinases qui sont liées entre autres choses à la conservation d'énergie sous forme d'ATP dans notre corps. Ces protéines, dont la famille Src fait partie, jouent un rôle dans beaucoup de fonctions cellulaire en permettant d'activer ou de désactiver des enzymes en leur attachant ou en leur détachant un ATP.

=== Maladies liées à SH2 ===

Les protéines du domaine SH2 ont été identifiées comme modifiées lors de maladies héréditaires ou sporadiques (sporadique: maladie  atteignant quelques individus isolés)

SH2 domain proteins have been identiﬁed as altered in hereditary or sporadic human diseases.  Known mutations in the genes for 18 distinct SH2 domain proteins contribute to human disorders, including cancers and leukemias, developmental disorders, diabetes, and immunodeﬁciencies.  These can arise from either loss- or gain-of-function mutations in SH2 domains.

Par exemple, dans le gène SH2D1A, sa mutation provoque un mauvais repliement de la protéine et un mauvais fonctionnement de celle-ci.

=== Structure 3D de SH2 ===

{{ :f20815:start:333083:sh2_domain.png?200 |}}

Sa structure est formée de deux α-hélices à l'extérieur et de sept ß-feuillets anti-parallèles à l'intérieur formant un large feuillet bêta. On peut le trouver dans de nombreux organismes, notamment l'être humain chez qui il est présent dans au moins 115 protéines identifiées à ce jour.

=== Interactions favorables entre SH2 et d'autres domaines ===

SH2 sont des éléments qui répondent habituellement à la phosphorylation sur tyrosine en liant la séquence phosphorylée. En tant que tels , ils constituent des éléments essentiels dans la régulation de la tyrosine kinase de processus cellulaires. En raison des interactions du domaine SH2 résultant de la phosphorylation, ces éléments fournissent un circuit réglable le long du quel des signaux peuvent être transmis dans les meilleurs délais.

==== SH3_1 Family ====
=== Introduction au domaine SH3 ===

Les domaines SH3 sont trouvés dans les protéines codant le génome humain. Son rôle principal est de réguler l'état d'activité des adaptateurs de protéines ainsi que des Tyrosine Kinases et on pense qu'il augmenterait la spécificité du substrat de certaines Tyrosine Kinases.  Les domaines SH3 classiques sont contraints chez les humains, aux protéines intracellulaires et donc à l'équilibre de Donnan (Cet équilibre sous-entend les règles d'échanges des ions au travers de la membrane plasmique. Les protéines intracellulaires sont chargées négativement (anions). La membrane cellulaire étant imperméable pour ces protéines, elles ne peuvent donc pas passer à l'extérieur de la cellule. Ces protéines sont, par conséquent, responsables de l' électronégativité de la cellule. Puisque elles sont électriquement négatives, elles retiennent les ions chargés positivement et essentiellement les ions potassium (cations)).

=== Maladies liées à SH3 ===

Des expressions réduites de SH3GL2 (SH3-domain GRB2-like 2, une protéine très présente dans le cerveau, la glande pituitaire et les reins) causées par des altérations moléculaires, sont impliquées dans la formation de tumeurs au niveau de la tête, le cou, les seins et les tissus gastriques. Mais ce sont plus souvent des cas de cancers sporadiques.

Par contre l'augmentation de l'expression du niveau de SH3GL2 dans les neurones est lié à l'activation de la kinase "stress" c-Jun qui provoque la mort du neurone. L'expression trop importante de SH3GL2 est maintenant considérée comme un nouvel indicateur de la progression de la maladie d'Alzheimer.


=== Structure 3D de SH3 ===
{{ :f20815:start:333083:sh3_domain.png?200 }}

Le repli SH3 consiste en 2 feuillets anti-parallèles bêta qui sont posés aux angles droits entre eux, composés de 5 ou 6 bêta . A l'intérieure du repli, il y a deux boucles variables appellées boucles RT et n-Src. Ce type de repli est un ancien repli retrouvé dans les eucaryotes et procaryotes.

=== Interactions favorables entre SH3 et d'autres domaines ===

Il est trouvé dans des protéines qui intéragissent et régulent un ensemble de protéines complexes. Il est susceptibles de se lier à des peptides riches en proline, se trouve en général à la surface des protéines cibles.

===== Création du profil PSSM=====
==== Alignement des séquences ====
Utilisation des sites CLUSTAL Omega et MUSCLE pour effectuer les alignements des différentes séquences.
On va sur les deux liens et on leur demande de nous aligner les fichiers "to-be-aligned" qu'on a récupéré lors de nos recherches.
====Parsage du fichier avec les séquences alignées====
On lit le fichier des séquences alignées qui contient un ">" avant chaque séquence. Il suffit de placer les séquences entre ">" dans une liste.

====Application de la Formule générale====

$m_{u,a} =  log(q_{u,a}/p_a)$

la formule de base de $q_{u,a} = (n_{u,a}+1)/(N_seq+20)$

mais en ajoutant les pseudocounts la formule devient$q_{u,a} = (alpha*f_{u,a}+beta*p_a)/(alpha+beta)$, elle permet de ne jamais obtenir un résultat = 0 et donc évite de se retrouver avec des log(0).

$p_a$ = au swissprot, un tableau de pourcentages donné dans les slides. (correspond à la probabilité qu'on trouve l'acide aminé à n'importe quelle position dans les séquences )

$alpha = N_seq - 1$

$beta = sqrt(N_seq)$

$N_{seq} = $ nombre de séquences total




===== Validation et présentation =====

==== Le weblogo ====
En allant sur le site proposé par les slides, les domaines choisis étaint considérés comme trop grands pour la création d'un weblogo, une recherche s'impose donc. C'est donc avec le site weblogo3 que nous allons effectuer le logo.

Pour pouvoir comparer le logo avec les résultats obtenus, le meilleur moyen et d'utiliser la séquence consensus. Cette séquence est obtenue en ne gardant que l'acide aminé avec le plus grand score pour chaque position.

=== Quelles sont les positions conservées ===
En comparant le weblogo avec les profils créés à priori toutes les positions sont conservées et correspondent bien au weblogo donné par le site.
Il se pourrait que sur des égalités de score les positions ne soient pas conservées mais ce, dans un très petit nombre de cas.

En comparant les deux PSSM et les deux weblogo, on peut se rendre compte des différences entre les algorithmes d'alignement multiples de MUSCLE et de Clustal. Alors que MUSCLE produit un alignement plus court et donc introduit moins de "gap", il semble moins efficace pour détecter les domaines.

=== Comparaison avec les informations de Pfam ===

En comparant le logo SH2 créé avec CLUSTAL, avec le logo HMM (http://pfam.xfam.org/family/PF00017#tabview=tab4) du profil Pfam on essaye de retrouver le domaine dans notre logo qu'on retrouve au niveau de la 2418 colonne et se termine en 2564, ce qui fait beaucoup de résidus pour un domaine d'une longueur de 76 qui là en fait 146 où une diminution du score se produit.

Si maintenant on compare le logo créé avec MUSCLE, on remarque que l'algorithme de MUSCLE produit un alignement plus court, de ce fait il introduit moins de "gap" mais semble moins efficace à détecter les domaines car celui-ci est tellement dispersé qu'il n'est pas identifiable.




====Bibliographie (simplifiée)====
  * http://pawsonlab.mshri.on.ca/index.php?option=com_content&task=view&id=178&Itemid=64
  * http://pawsonlab.mshri.on.ca/index.php?option=com_content&task=view&Itemid=64&id=179
  * http://pfam.xfam.org/family/SH2
  * http://pfam.xfam.org/clan/SH3
  * http://www.univ-montp1.fr/recherche/unites_de_recherche/physiologie_medecine_experimentale_du_caeur_et_des_muscles_inserm_u1046/equipes_teams/pour_en_savoir_plus_more_details/repertoire_de_a_a_z_de_divers_domaines_au_sein_des_proteines/le_domaine_sh2
  * http://www.univ-montp1.fr/recherche/unites_de_recherche/physiologie_medecine_experimentale_du_caeur_et_des_muscles_inserm_u1046/equipes_teams/pour_en_savoir_plus_more_details/repertoire_de_a_a_z_de_divers_domaines_au_sein_des_proteines/le_domaine_sh3
  * http://en.wikipedia.org/wiki/SH2_domain
  * http://en.wikipedia.org/wiki/SH3_domain
  * http://pfam.xfam.org/family/PF00017
  * http://pfam.xfam.org/family/PF00018
  * http://es.wikipedia.org/wiki/Dominio_SH3
  * http://www.ncbi.nlm.nih.gov/pubmed/7542925
  * http://atlasgeneticsoncology.org/Genes/SH3GL2ID44345ch9p22.html
  * http://lyrobossite.free.fr/Depolarisation_I_Le%20potentiel%20de%20repos_1.htm
