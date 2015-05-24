=======Introduction à la bioinformatique: portfolio=======


======Projet 1: Implémentation de l’algorithme Needleman-Wunsch et l’algorithme de Smith Waterman qui sera comparée avec le logiciel LALIGN.======

{{:f20815:start:333083:mini_projet_1.zip|}}


====Partie 1.1 Implémentation de deux ADT (Sequence, Score)====
Le but de ces des classes est de pouvoir extraire d'un fichier soit des sequences ou des matrices de substitution et de pouvoir les manipuler.

Ce qui va nous intéresser pour les séquences va être de bien parser le fichier et d'ignorer toute l'information superflue.

Et pour les matrices de substitution, en stockant celles-ci avec leur liste d'acides aminés dans l'ordre on va pouvoir renvoyer le score pour une paire d'acides.

====Partie 1.2 Implémentation de l'algorithme de Needleman-Wunsch (Alignement global avec pénalité affine)====

On commencera par parler de l'algorithme, ses méthodes et on finira avec une comparaison entre mes résultats et ceux de LALIGN.
Il y a 4 parties intéressantes à l'algorithme: la pénalité affine, l'initialisation des trois matrices et leur utilité, les remplissage de ces 3 matrices et la fonction récursive qui effectue un backtraking du chemin effectué.

===La pénalité affine===
Le but de la pénalité affine est de calculer le score de gap, pour se faire on doit pouvoir calculer la probabilité qu'un gap ait été inséré. Les deux matrices supplémentaires ajoutées vont permettre de calculer cette probabilité supplémentaire.

S = matrice des scores d'alignement,
V = matrice des scores de gap dans la 1ère séquence,
W = matrice des scores de gap dans la 2ème séquence,
I = pénalité de début d'un gap,
E = pénalité d'extension d'un gap,
T = matrice de substitution.

$$ V_{i, j} = max(V_{i-1, j} - E, S_{i-1, j} - I - E) $$


$$ S_{i, j} = max(S_{i-1, j-1 + T_{i, j}, V_{i, j}, W_{i, j}) $$

===L'initialisation des trois matrices===
On va initialiser les 3 matrices avec -inf car c'est la seule valeur qui va nous permettre d'avoir des bonnes comparaisons lorsqu'on recherchera le score maximum. Les matrices sont plus longues d'une case pour permettre de calculer le score en début de séquence à partir de valeurs par défaut qui sont déterminées à partir des pénalités. Ces valeurs par défaut permettent un gap en début de séquence.


<code python>
self.S = [i[:] for i in [[float("-inf")]*(len(sequence2)+1)]*(len(sequence1)+1)]
self.S[0][0] = 0

self.V = [i[:] for i in [[float("-inf")]*(len(sequence2)+1)]*(len(sequence1)+1)]

self.W = [i[:] for i in [[float("-inf")]*(len(sequence2)+1)]*(len(sequence1)+1)]

for i in range(1,len(self.sequence1)+1):
	self.S[i][0] = -gapStart-(i-1)*gapExtend
	self.V[i][0] = self.S[i][0]
	self.W[i][0] = self.S[i][0]

for j in range(1,len(self.sequence2)+1):
	self.S[0][j] = -gapStart-(j-1)*gapExtend
	self.V[0][j] = self.S[0][j]
	self.W[0][j] = self.S[0][j]

</code>

===Remplissage des trois matrices===

<code python>
def fillMatrix(self):
	for i in range(1,len(self.sequence1)+1):
		for j in range(1,len(self.sequence2)+1):
			self.V[i][j] = max(self.S[i-1][j]-self.gapStart-self.gapExtend, self.V[i-1][j]-self.gapExtend)
			self.W[i][j] = max(self.S[i][j-1]-self.gapStart-self.gapExtend, self.W[i][j-1]-self.gapExtend)
			self.S[i][j] = max(self.S[i-1][j-1] + self.score[self.sequence1[i-1],self.sequence2[j-1]], self.V[i][j], self.W[i][j]) # ATTENTION AVEC SCORE !! [][] ou [A,B]
			if self.alignmentType == "local":
				self.S[i][j] = max(self.S[i][j],0)
				if self.S[i][j] > self.S[self.lengths[0]][self.lengths[1]]:
					self.lengths = [i,j]
</code>

===Backtracking===
On va parcourir en sens inverse les chemins partant du score global optimal. La condition d'arrêt est l'arrivée à la case [0][0] pour l'alignement global et n'importe quelle case valant 0 pour l'alignement local. On crée le chemin en suivant les scores de gap et d'alignement pour savoir dans quelle direction se déplacer.

<code python>
def traceback(self,i,j, alignment1,alignment2):
	if ((self.alignmentType == "global" and (i > 0 or j > 0)) or (self.alignmentType == "local" and self.S[i][j] > 0)):
		if (self.isFromDiagonal(i, j)):
			self.traceback(i-1, j-1, self.sequence1[i-1] + alignment1, self.sequence2[j-1] + alignment2)
		if (self.isFromLeft(i, j)):
			self.traceback(i-1, j, self.sequence1[i-1] + alignment1, "-" + alignment2)
		if (self.isFromUp(i, j)):
			self.traceback(i, j-1, "-" + alignment1, self.sequence2[j-1] + alignment2)
		else:
			self.finalAlignment.append((alignment1, alignment2))
</code>

===Comparaison des résultats obtenus avec LALIGN===

Je testerai avec les 2 premières séquences du fichier SH2-sequences.fasta

==Alignement global==

- Mes résultats:

<code>
Welcome to the sequence aligner,
 Introduce the following information to obtain the alignment between 2 sequences.
Enter the 1st sequence filename: data MP1 2015/SH2-sequences.fasta
Choose 1 sequence from all given in the file ( from 0 to 5 ): 0
Same thing for the 2nd sequence,
Enter the 2nd sequence filename: data MP1 2015/SH2-sequences.fasta
Choose 1 sequence from all given in the file ( from 0 to 5 ): 1
Enter the filename of the score matrix you want to use: ScoreMatrix/blosum62.txt
Enter the opening gap penalty: 5
Enter the extending gap penalty: 2
Give the alignmentType, 'global' or 'local': global
YFGKLGRKDAERQLLSFGNPRGTFLIRESETTKGAYSLSIRDWDDMKGDHVKHYKIRKLDNGGYYITTRAQFETLQQLVQHYSERAAGLCCRLVVPC
::::..:...::.::...:::::::.::::::::::.::..:.:..::..::::::::::.::.:::.:.::..:::::..::..:.:::.::...:
YFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHYKIRKLDSGGFYITSRTQFNSLQQLVAYYSKHADGLCHRLTTVC
66.0% identity
100.0% similarity
0.0% gap
Length : 97
score : 353
</code>

-Résultats de LALIGN:

<code>
The best scores are:                                      n-w bits E(1)
seq2 98 bp                                         (  98)  353 66.8  7e-179

>>seq2 98 bp                                              (98 aa)
 n-w opt: 353  Z-score: 334.9  bits: 66.8 E(1): 7e-179
global/global (N-W) score: 353; 66.3% identity (87.8% similar) in 98 aa overlap (1-98:1-98)

               10        20        30        40        50        60
seq1   WYFGKLGRKDAERQLLSFGNPRGTFLIRESETTKGAYSLSIRDWDDMKGDHVKHYKIRKL
       :::::. :...:: ::.  :::::::.:::::::::: ::. :.:. :: .:::::::::
seq2   WYFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHYKIRKL
               10        20        30        40        50        60

               70        80        90        
seq1   DNGGYYITTRAQFETLQQLVQHYSERAAGLCCRLVVPC
       :.::.:::.:.::..::::: .::..: ::: ::.. :
seq2   DSGGFYITSRTQFNSLQQLVAYYSKHADGLCHRLTTVC
               70        80        90
</code>

==Alignement local==
-Mes résultats
<code>
Welcome to the sequence aligner,
 Introduce the following information to obtain the alignment between 2 sequences.
Enter the 1st sequence filename: data MP1 2015/SH2-sequences.fasta
Choose 1 sequence from all given in the file ( from 0 to 5 ): 0
Same thing for the 2nd sequence,
Enter the 2nd sequence filename: data MP1 2015/SH2-sequences.fasta
Choose 1 sequence from all given in the file ( from 0 to 5 ): 1
Enter the filename of the score matrix you want to use: ScoreMatrix/blosum62.txt
Enter the opening gap penalty: 5
Enter the extending gap penalty: 1
Give the alignmentType, 'global' or 'local': local
YFGKLGRKDAERQLLSFGNPRGTFLIRESETTKGAYSLSIRDWDDMKGDHVKHYKIRKLDNGGYYITTRAQFETLQQLVQHYSERAAGLCCRLVVPC
::::..:...::.::...:::::::.::::::::::.::..:.:..::..::::::::::.::.:::.:.::..:::::..::..:.:::.::...:
YFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHYKIRKLDSGGFYITSRTQFNSLQQLVAYYSKHADGLCHRLTTVC
66.0% identity
100.0% similarity
0.0% gap
Length : 97
score : 353
</code>


-Résultats de LALIGN:

<code>
Waterman-Eggert score: 353;  85.1 bits; E(1) <  2.3e-22
66.3% identity (87.8% similar) in 98 aa overlap (1-98:1-98)

               10        20        30        40        50        60
seq1   WYFGKLGRKDAERQLLSFGNPRGTFLIRESETTKGAYSLSIRDWDDMKGDHVKHYKIRKL
       :::::. :...:: ::.  :::::::.:::::::::: ::. :.:. :: .:::::::::
seq2   WYFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHYKIRKL
               10        20        30        40        50        60

               70        80        90        
seq1   DNGGYYITTRAQFETLQQLVQHYSERAAGLCCRLVVPC
       :.::.:::.:.::..::::: .::..: ::: ::.. :
seq2   DSGGFYITSRTQFNSLQQLVAYYSKHADGLCHRLTTVC
               70        80        90        
</code>
