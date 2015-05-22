=====Projet 3: La construction de profils=====

==== Introduction ====


Une PSSM (Position Specific Scoring Matrix) permet de représenter un motif dans un groupe de séquences (de protéines dans notre cas). Elles sont très utiles pour représenter une famille ou un domaine de protéines dans son ensemble.
Ce troisième projet consiste en l'implémentation d'un algorithme de création d'une PSSM et dans la vérification de sa pertinence biologique. Le principe est semblable à la construction d'une matrice de substitution mais en calculant le score pour chaque position de chaque acide aminé dans le groupe de séquence donné. On travaille avec un alignement multiple de séquences provenant du même domaine ou de la même famille dans le but de créer une PSSM spécifique à ce domaine ou cette famille.

Pour commencer, je vais faire une introduction des deux domaines utilisés, s'en suivra une présentation de mon algorithme et de mon implémentation.

L'archive zip est disponible ici :

=== Domaine SH2 ===

Le premier domaine que j'ai analysé est le domaine SH2. Ce domaine est contenu principalement dans l'oncoprotéine Src dont il tire son nom Src Homologuous 2 mais également dans d'autre protéines chargées de l'échange de messages entre les cellules.

Le domaine SH2 permet à une protéine de s'attacher à d'autre protéines par un résidu tyrosine, c'est pourquoi on le trouve souvent dans des protéines impliquées dans la transduction d'un signal par les récepteurs tyrosines kinases qui sont liées entre autres choses à la conservation d'énergie sous forme d'ATP dans notre corps. Ces protéines, dont la famille Src fait partie, jouent un rôle dans beaucoup de fonctions cellulaire en permettant d'activer ou de désactiver des enzymes en leur attachant ou en leur détachant un ATP.
On le trouve également dans des protéines impliquées dans la transcription ou dans des enzymes commet les phosphatases dans lesquelles son but est de détecter des tyrosines.

Sa structure est formée de deux α-hélices à l'extérieur et de sept ß-feuillets anti-parallèles à l'intérieur. On peut le trouver dans de nombreux organismes, notamment l'être humain chez qui il est présent dans au moins 115 protéines identifiées à ce jour.


Petite introduction sur les domaines avec lesquels on va travailler dans ce projet.
Le domaine SH2 
Dans SH2, es intéraction protéine-protéine jouent un rôle important dans le développement et la croissance cellulaire. Les domaines modulaires, sous-unités des protéines, régulent ces intéractions en identifiant de courtes séquences peptidiques.
Un des domaines les plus proéminents est le domaine SH2 (100 acides aminés, trouvé dans 111 protéines hmaines) qui joue un rôle très important dans la communication cellulaire.

Le SH3 est un domaine moins important (60 acides aminés et 300 domaines SH3 sont trouvés dans les protéines codant le génome humain). Son rôle principal est de réguler l'état d'activité des adaptateurs de protéines ainsi que des Tyrosine Kinases et on pense qu'il augmenterait la spécificité du substrat de certaines Tyrosine Kinases. Il est trouvé dans des protéines qui intéragissent et régulent un ensemble de protéines complexes. Il est susceptibles de se lier à des peptides riches en proline, se trouve en général à la surface des protéines cibles. Les domaines SH3 classiques sont contraints chez les humains, aux protéines intracellulaires.

== Relations interdomaines ==
Selon le site Pfam, SH2 aurait 8 intéractions:
Pkinase_Tyr, Y_phosphatase, SOCS_box, SH2, SH3_1, SH3_2, STAT_bind et C1_1 
Et SH3 aurait 16 interactions:
SH2, SH3_1, Ank, PH, ubiquitin, LIM, B1,RhoGEF, P53, Pkinase_Tyr, SH2, SH3_1, Guanylate_kin, SH3_2, SH3_5

== Implications dans certaines maladies ? ==
SH2: CSK in Human Diseases

Csk's interaction with a phosphotase ("Lyp", gene product of PTPN22) is possibly associated with the increased autoimmune diseases associated with PTPN22 mutations.

SH3: Cdc25s in Human Disease

The Cdc25s, and in particular Cdc25A and Cdc25B, are proto-oncogenes in humans and have been shown to be overexpressed in a number of cancers.The central role of Cdc25s in the cell cycle has garnered them considerable attention from the pharmaceutical industry as potential targets for novel chemotherapeutic. To date, no clinically viable compounds targeting these enzymes have been described.


== Structure 3D de ces proteines ==
SH2: La structure consiste en un large feuillet bêta(composé de 7 bêta-brins) pris de revers par 2 hélices alpha.

SH3: Le repli SH3 consiste en 2 feuillets ant-parallèles bêta qui sont posés aux angles droits entre eux. A l'intérieure du repli, il y a deux boucles variables appellées boucles RT et n-Src. Ce type de repli est un ancien repli retrouvé dans les eucaryotes et procaryotes.

=== Travail de recherche ===
Pour trouver les séquences de ces domaines, j'ai utilisé le site d'Uniprot principalement. J'ai voulu utiliser PFam mais j'ai eu plus de mal à trouver mon bonheur.

  - Recherche effectuée: domain:SH2 AND organism:human 
  - Recherche effectuée: domain:SH3 AND organism:human AND gene name:SH3.


==== Création du profil PSSM====
=== Alignement des séquences ===
Utilisation des sites CLUSTAL Omega et MUSCLE pour effectuer les alignements des différentes séquences.
On va sur les deux liens et on leur demande de nous aligner les fichiers "to-be-aligned" qu'on a récupéré lors de nos recherches.
===Parsage du fichier avec les séquences alignées===
On lit le fichier des séquences alignées qui contient un ">" avant chaque séquence. Il suffit de placer les séquences entre ">" dans une liste.

===Application de la Formule générale===

$m_{u,a} =  log(q_{u,a}/p_a)$

la formule de base de $q_{u,a} = (n_{u,a}+1)/(N_seq+20)

mais en ajoutant les pseudocounts la formule devient$q_{u,a} = (alpha*f_{u,a}+beta*p_a)/(alpha+beta)$, elle permet de ne jamais obtenir un résultat = 0 et donc évite de se retrouver avec des log(0).

$p_a$ = au swissprot, un tableau de pourcentages donné dans les slides. (correspond à la probabilité qu'on trouve l'acide aminé à n'importe quelle position dans les séquences )

$alpha = N_seq - 1$

$beta = sqrt(N_seq)$

$N_{seq} = nombre de séquences total$




==== Validation et présentation ====

=== Le weblogo ===
En allant sur le site proposé par les slides, les domaines choisis étaint considérés comme trop grands pour la création d'un weblogo, une recherche s'impose donc. C'est donc avec le site weblogo3 que nous allons effectuer le logo.

Pour pouvoir comparer le logo avec les résultats obtenus, le meilleur moyen et d'utiliser la séquence consensus. Cette séquence est obtenue en ne gardant que l'acide aminé avec le plus grand score pour chaque position.

== Quelles sont les positions conservées ==
En comparant le weblogo avec les profils créés à priori toutes les positions sont conservées et correspondent bien au weblogo donné par le site.
Il se pourrait que sur des égalités de score les positions ne soient pas conservées mais ce dans un très petit nombre de cas.
== Est-ce que les acides aminés qui font partie des structures secondaires sont bien conservés ==

== Comparaison entre nos résultats et l'information obtenue sur le site PFAM pour les deux domaines ==

Je comprends pas comment on est sensés comparer de manière efficace ?

====Bibliographie (simplifiée)====
  * http://pawsonlab.mshri.on.ca/index.php?option=com_content&task=view&id=178&Itemid=64
  * http://pawsonlab.mshri.on.ca/index.php?option=com_content&task=view&Itemid=64&id=179
  * http://pfam.xfam.org/family/SH2
  * http://pfam.xfam.org/clan/SH3
  * http://www.univ-montp1.fr/recherche/unites_de_recherche/physiologie_medecine_experimentale_du_caeur_et_des_muscles_inserm_u1046/equipes_teams/pour_en_savoir_plus_more_details/repertoire_de_a_a_z_de_divers_domaines_au_sein_des_proteines/le_domaine_sh2
  * http://www.univ-montp1.fr/recherche/unites_de_recherche/physiologie_medecine_experimentale_du_caeur_et_des_muscles_inserm_u1046/equipes_teams/pour_en_savoir_plus_more_details/repertoire_de_a_a_z_de_divers_domaines_au_sein_des_proteines/le_domaine_sh3
  * http://en.wikipedia.org/wiki/SH2_domain
  * http://en.wikipedia.org/wiki/SH3_domain


