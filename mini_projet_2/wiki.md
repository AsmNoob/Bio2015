====== Projet 2: Une implémentation de l’algorithme BLOSUM construite en utilisant des données concernant les domaines protéiques BRD et PTB. ======

{{:f20814:start:333083:code.zip|mini_projet_2}}

==== 1. Exemple du fonctionnement de la méthode BLOSUM en utilisant une des familles ====
== 1.1 Si vous ne possédez pas les fichiers fasta vous devez d'abord les télécharger(ci-dessous). Renommer les fichiers en PR00109.txt et PR00401.txt. ==

 {{f20814:start:333083:fichiers.zip|fichiers}} 

== 1.2 Après avoir téléchargé les fichiers si ce n'était pas déjà fait, lancer la commande suivante dans le terminal: "python3 Parseur.py". Ceci va créer des fichiers contenant tous les différents BLOCKS de chaque famille. ==
== 1.3 Vous pouvez maintenant lancer la commande "python3 BLOSUM.py" ==
== 1.4 Explication de chaque étape du processus ==

1) On parcourt l'ensemble des blocks qui constituent la famille (A->E)
2) Pour chacun des blocs de séquences, on cherche la liste des groupes que l'ont peut trouver.

<code python>
listeGroupes = trouverGroupesDifferents(listeSequences,pourcentage)

#Place les sequences dans des groupes différents selon leur similitudes
def trouverGroupesDifferents(listeSequences, pourcentage):
	listeGroupes = []
	for sequence in listeSequences:
		if(listeGroupes == []):
			liste = [sequence]
			listeGroupes.append(liste)
		else:
			groupeTrouve = False
			for groupe in listeGroupes:
				if(not groupeTrouve):
					for i in range(len(groupe)):
						if(comparaisonSequences(sequence,groupe[i], pourcentage) and not groupeTrouve):
							groupe.append(sequence)
							groupeTrouve = True
			if(not groupeTrouve):
				listeM = [sequence]
				listeGroupes.append(listeM)
				groupeTrouve = True
	return listeGroupes

#Compare 2 sequences et si leur taux d'dentité est supérieur au pourcentage donné.
def comparaisonSequences(seq1,seq2,pourcentage):
	elementDeSimilitude = 0
	longueurSequence = len(seq1)
	res = False
	for i in range(longueurSequence):
		if(seq1[i] == seq2[i]):
			elementDeSimilitude+=1
	if((elementDeSimilitude/longueurSequence) >= pourcentage):
		res = True
	return res

</code>

3) Après avoir trouvé les différents groupes constituant le bloc, on cherche la matrice des fréquences pondérées qui correspond.

<code python>
matriceFrequencesPonderees(listeGroupes,['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])

#Calcul de la matrice des frequences pondérées
def matriceFrequencesPonderees(listeGroupes,listeProteines):
	matrice = []
	for i in range(len(listeProteines)):
		matrice.append([0 for m in range(len(listeProteines))])
		for j in range(len(listeProteines)):
			matrice[i][j] = calculFrequence(listeProteines[i],listeProteines[j],listeGroupes)
	return matrice

def colonneProteines(listeGroupes,indice):
	colonne = []
	for groupe in listeGroupes:
		for sequence in groupe:
			colonne.append(sequence[indice])
	return colonne

#Calcule la fréquence pondérée pour 2 proteines données 
def calculFrequence(proteine1,proteine2,listeGroupes):
	somme = 0
	for indice in range(len(listeGroupes[0][0])): #taille d'une sequence
		colonne = colonneProteines(listeGroupes,indice)#colonne de proteines sur tous les groupes
		for i in range(len(colonne)):
			if(colonne[i] == proteine1):
				if(proteine1 == proteine2):
					for j in range(i,len(colonne)):
						if(j != i and colonne[j] == proteine2 and calculGroupe(i,listeGroupes) != calculGroupe(j,listeGroupes)):
							somme+=(calculPoidsGroupe(i,listeGroupes)*(calculPoidsGroupe(j,listeGroupes)))
				else:
					for j in range(len(colonne)):
						if(j != i and colonne[j] == proteine2 and calculGroupe(i,listeGroupes) != calculGroupe(j,listeGroupes)):
							somme+=(calculPoidsGroupe(i,listeGroupes)*(calculPoidsGroupe(j,listeGroupes)))
	return somme

def calculGroupe(indiceColonne,listeGroupes):
	numeroGroupe = 0
	positionActuelle = len(listeGroupes[numeroGroupe]) #taille du premier groupe
	for groupe in listeGroupes:
		if(indiceColonne >= positionActuelle):
			numeroGroupe+=1
			positionActuelle+= len(listeGroupes[numeroGroupe])
	return numeroGroupe

#calcul du poids de la sequence, en utilisant son indice ds la colonne et la taille des différents groupes
def calculPoidsGroupe(indiceColonne,listeGroupes):
	numeroGroupe = 0
	positionActuelle = len(listeGroupes[numeroGroupe]) #taille du premier groupe
	for groupe in listeGroupes:
		if(indiceColonne >= positionActuelle):
			numeroGroupe+=1
			positionActuelle+= len(listeGroupes[numeroGroupe])
	return (1/(len(listeGroupes[numeroGroupe])))


</code>

4) On place ensuite la matrice trouvée dans une liste. On effectue les trois points précédents pour tous les blocs de la famille et on effectue une somme normalisée de l'ensemble des matrices trouvées

<code python>
matriceResultante = moyenneMatrices(listeMatrices)

#Moyenne des matrices passées en paramètres
def moyenneMatrices(listeMatrices):
	matrice = [[0 for i in range(20)] for j in range(20)]
	for i in range(len(listeMatrices[0])):
		for j in range(len(listeMatrices[0][0])):
			for n in range(len(listeMatrices)):

				matrice[i][j] += listeMatrices[n][i][j]
			matrice[i][j] = matrice[i][j]/len(listeMatrices)
	return matrice

</code>

5) On calcule ensuite la matrice des probabilités d'occurences:

<code python>
matriceOccurences = calculProbabiliteOccurence(matriceResultante)

def calculProbabiliteOccurence(matrice):
	sommePartieSuperieure = calculPartieSuperieureMatrice(matrice)
	if(sommePartieSuperieure != 0):
		for i in range(len(matrice)):
			for j in range(len(matrice[0])):
				matrice[i][j]= (matrice[i][j]/sommePartieSuperieure)
	return matrice

def calculPartieSuperieureMatrice(matrice):
	somme = 0
	for i in range(len(matrice)):
		for j in range(len(matrice[0])):
			if(j >= i):
				somme+=matrice[i][j]
	return somme
</code>

6) Et on finit avec la matrice finale:

<code python>
matriceFinale = calculMatriceFinale(matriceOccurences)

def calculMatriceFinale(matriceOccurences):
	copieMatrice = copy.deepcopy(matriceOccurences)
	for i in range(len(matriceOccurences)):
		for j in range(len(matriceOccurences[0])):
			matriceOccurences[i][j] = calculTauxLogChance(copieMatrice,i,j)
	return matriceOccurences


def calculTauxLogChance(matriceOccurences,i,j):
	try:
		if((matriceOccurences[i][j])/(frequencePrevuePourAlignement(matriceOccurences,i,j)) == 0):
			res = 0
		else:
			res = 2*(log((matriceOccurences[i][j])/(frequencePrevuePourAlignement(matriceOccurences,i,j)), 2))
	except(ZeroDivisionError):
		res = 0
	return res

def frequencePrevuePourAlignement(matrice,indiceProteine1,indiceProteine2):
	if(indiceProteine1 == indiceProteine2):
		frequence = frequencePrevueParResidu(matrice,indiceProteine2)**2
	else:
		frequence = (2*frequencePrevueParResidu(matrice,indiceProteine1)*frequencePrevueParResidu(matrice,indiceProteine2))
	return frequence

def frequencePrevueParResidu(matriceOccurences,indiceProteine):
	somme = 0
	for i in range(len(matriceOccurences)):
		if(i != indiceProteine):
			somme+=matriceOccurences[i][indiceProteine]
	somme= somme/2
	somme+=matriceOccurences[indiceProteine][indiceProteine]
	return somme
	

</code>

==== 2. Similarités entre les matrices trouvées et la matrice BLOSUM62 ====
Tout d'abord un affichage des 4 matrices trouvées:
<code python>
Substitution matrix generated with blocks from the Tyrosine kinase family
Identity% = 40.0%
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y   
A   -4  -2  -4  1   -8  3   0   -3  0   4   -2  -2  -2  -1  1   2   0   -2  0   -4  
C   -2  7   0   -10 0   0   0   -10 0   0   0   0   2   0   -13 -10 -2  -9  -11 0   
D   -4  0   3   3   0   3   -3  -11 -3  -2  0   0   -6  2   2   -2  -5  -7  0   -11 
E   1   -10 3   3   -11 1   1   -5  1   1   -1  2   -6  4   -2  2   -4  -5  -12 -7  
F   -8  0   0   -11 9   -4  0   -2  0   -7  -6  -2  0   0   -7  -5  -10 -6  0   -4  
G   3   0   3   1   -4  0   2   -8  5   -2  1   1   -7  4   -1  -1  -2  -3  -3  0   
H   0   0   -3  1   0   2   0   0   0   -4  0   6   -6  0   4   0   -3  0   0   0   
I   -3  -10 -11 -5  -2  -8  0   2   -9  3   8   0   -11 0   -10 -10 -6  3   0   0   
K   0   0   -3  1   0   5   0   -9  0   -1  0   0   -8  0   2   4   3   -7  0   0   
L   4   0   -2  1   -7  -2  -4  3   -1  -7  -5  -3  -6  -3  -5  3   1   5   0   5   
M   -2  0   0   -1  -6  1   0   8   0   -5  0   0   0   0   -7  0   -3  1   0   0   
N   -2  0   0   2   -2  1   6   0   0   -3  0   0   0   -6  0   -4  -1  0   0   6   
P   -2  2   -6  -6  0   -7  -6  -11 -8  -6  0   0   5   -9  -5  -2  -2  -8  0   0   
Q   -1  0   2   4   0   4   0   0   0   -3  0   -6  -9  0   2   -1  -2  0   0   0   
R   1   -13 2   -2  -7  -1  4   -10 2   -5  -7  0   -5  2   4   -1  -1  0   0   1   
S   2   -10 -2  2   -5  -1  0   -10 4   3   0   -4  -2  -1  -1  -8  3   -5  -4  -1  
T   0   -2  -5  -4  -10 -2  -3  -6  3   1   -3  -1  -2  -2  -1  3   6   -2  0   -5  
V   -2  -9  -7  -5  -6  -3  0   3   -7  5   1   0   -8  0   0   -5  -2  5   0   5   
W   0   -11 0   -12 0   -3  0   0   0   0   0   0   0   0   0   -4  0   0   9   -9  
Y   -4  0   -11 -7  -4  0   0   0   0   5   0   6   0   0   1   -1  -5  5   -9  -5  
Substitution matrix generated with blocks from the Tyrosine kinase family
Identity% = 70.0%
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y   
A   2   0   0   1   -5  -1  2   -3  1   -3  -2  1   0   1   -1  1   1   -2  -7  -1  
C   0   7   -6  -8  -13 -9  -8  -4  -7  -5  -6  -6  -5  -6  -8  -4  0   -5  -12 -10 
D   0   -6  4   2   -8  0   -1  -4  2   -3  -3  3   -3  1   -1  1   -1  -3  -12 -3  
E   1   -8  2   4   -4  1   1   -4  2   -3  -3  2   -3  2   -1  1   -1  -3  -6  -2  
F   -5  -13 -8  -4  8   -7  -3  -1  -5  -2  -3  -7  -10 -6  -7  -7  -5  -1  -7  2   
G   -1  -9  0   1   -7  9   -1  -6  0   -4  -5  0   -6  0   -2  0   0   -5  -7  -4  
H   2   -8  -1  1   -3  -1  4   -5  2   -2  -5  2   -4  3   0   0   0   -4  -5  3   
I   -3  -4  -4  -4  -1  -6  -5  5   -3  4   0   -4  -3  -3  -4  -5  -4  5   -6  -3  
K   1   -7  2   2   -5  0   2   -3  2   -1  -3  2   -5  2   1   1   1   -2  -6  -2  
L   -3  -5  -3  -3  -2  -4  -2  4   -1  5   0   -1  -5  -1  -2  -3  -2  4   -7  -2  
M   -2  -6  -3  -3  -3  -5  -5  0   -3  0   7   -3  -11 -1  -4  -3  -3  0   -8  -5  
N   1   -6  3   2   -7  0   2   -4  2   -1  -3  3   -5  2   0   1   0   -3  -8  -3  
P   0   -5  -3  -3  -10 -6  -4  -3  -5  -5  -11 -5  6   -4  -5  1   -2  -5  -11 -7  
Q   1   -6  1   2   -6  0   3   -3  2   -1  -1  2   -4  3   1   1   -1  -3  -6  -2  
R   -1  -8  -1  -1  -7  -2  0   -4  1   -2  -4  0   -5  1   5   -1  -2  -4  -9  -2  
S   1   -4  1   1   -7  0   0   -5  1   -3  -3  1   1   1   -1  3   3   -3  -9  -4  
T   1   0   -1  -1  -5  0   0   -4  1   -2  -3  0   -2  -1  -2  3   6   0   -9  -2  
V   -2  -5  -3  -3  -1  -5  -4  5   -2  4   0   -3  -5  -3  -4  -3  0   6   -5  -3  
W   -7  -12 -12 -6  -7  -7  -5  -6  -6  -7  -8  -8  -11 -6  -9  -9  -9  -5  9   -4  
Y   -1  -10 -3  -2  2   -4  3   -3  -2  -2  -5  -3  -7  -2  -2  -4  -2  -3  -4  7   
Substitution matrix generated with blocks from the SH2 family
Identity% = 40.0%
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y   
A   5   -6  4   1   -9  -2  -7  -10 -5  2   -1  -2  -1  0   -5  -1  4   -4  -8  0   
C   -6  0   -3  0   4   -7  -9  -8  -5  0   -3  0   -7  1   3   -1  -8  -1  -4  5   
D   4   -3  3   2   -11 3   -6  -3  -1  -1  -4  1   3   0   -5  1   2   -6  -9  -8  
E   1   0   2   5   -6  2   -6  -4  -1  -3  -3  0   -4  -1  1   1   -1  -1  3   -4  
F   -9  4   -11 -6  6   -8  0   -10 -11 -8  -7  0   0   -11 -10 -9  -18 -3  -3  5   
G   -2  -7  3   2   -8  7   -10 -9  -5  -7  0   2   2   -2  0   2   -1  -10 -9  -7  
H   -7  -9  -6  -6  0   -10 7   -7  -6  -8  -11 -3  -8  -6  -7  -1  -2  -9  -9  -4  
I   -10 -8  -3  -4  -10 -9  -7  4   1   1   2   -1  -3  -1  -2  -5  -1  3   -8  -9  
K   -5  -5  -1  -1  -11 -5  -6  1   3   0   3   0   0   3   2   -1  2   -4  -2  -7  
L   2   0   -1  -3  -8  -7  -8  1   0   4   3   -3  -9  1   -1  -4  2   -3  -7  1   
M   -1  -3  -4  -3  -7  0   -11 2   3   3   3   0   1   0   0   -5  4   -6  -7  -5  
N   -2  0   1   0   0   2   -3  -1  0   -3  0   7   -1  1   0   0   -3  -3  4   -2  
P   -1  -7  3   -4  0   2   -8  -3  0   -9  1   -1  5   -1  -8  0   -5  2   -2  -8  
Q   0   1   0   -1  -11 -2  -6  -1  3   1   0   1   -1  2   3   -4  1   -5  -7  -1  
R   -5  3   -5  1   -10 0   -7  -2  2   -1  0   0   -8  3   3   1   0   -3  2   0   
S   -1  -1  1   1   -9  2   -1  -5  -1  -4  -5  0   0   -4  1   4   -1  -1  1   -8  
T   4   -8  2   -1  -18 -1  -2  -1  2   2   4   -3  -5  1   0   -1  3   -4  -5  -11 
V   -4  -1  -6  -1  -3  -10 -9  3   -4  -3  -6  -3  2   -5  -3  -1  -4  4   -1  -4  
W   -8  -4  -9  3   -3  -9  -9  -8  -2  -7  -7  4   -2  -7  2   1   -5  -1  11  -5  
Y   0   5   -8  -4  5   -7  -4  -9  -7  1   -5  -2  -8  -1  0   -8  -11 -4  -5  5   
Substitution matrix generated with blocks from the SH2 family
Identity% = 70.0%
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y   
A   6   -2  1   0   -6  -1  -2  -4  0   -5  1   0   1   0   -1  0   0   -1  -5  -7  
C   -2  3   1   -3  1   -3  -2  0   0   -1  0   0   -1  -1  0   3   0   -2  -3  2   
D   1   1   3   2   -5  0   -1  -3  1   -5  -1  2   1   2   -1  0   0   -3  -4  -6  
E   0   -3  2   5   -9  -1  -3  -2  1   -5  -1  0   -1  2   -2  -1  -1  -3  -6  -7  
F   -6  1   -5  -9  7   -5  0   -7  -4  -3  -2  -4  -7  -7  -4  -7  -6  -2  0   3   
G   -1  -3  0   -1  -5  6   -5  -4  0   -6  -5  1   0   0   0   0   -1  -6  -5  -5  
H   -2  -2  -1  -3  0   -5  6   -6  -1  -5  -1  0   -2  -1  -3  -2  -1  -3  -1  3   
I   -4  0   -3  -2  -7  -4  -6  6   -2  2   -1  -3  0   -2  -2  -4  -3  3   -7  -5  
K   0   0   1   1   -4  0   -1  -2  2   -3  0   1   1   2   1   0   1   -3  -4  -4  
L   -5  -1  -5  -5  -3  -6  -5  2   -3  5   3   -5  -2  -3  -2  -5  -5  -2  -4  -4  
M   1   0   -1  -1  -2  -5  -1  -1  0   3   2   0   0   1   1   -2  -1  -2  -2  -4  
N   0   0   2   0   -4  1   0   -3  1   -5  0   2   2   2   0   1   1   -2  -3  -4  
P   1   -1  1   -1  -7  0   -2  0   1   -2  0   2   5   2   -3  1   0   -1  -3  -5  
Q   0   -1  2   2   -7  0   -1  -2  2   -3  1   2   2   2   0   0   1   -2  -3  -4  
R   -1  0   -1  -2  -4  0   -3  -2  1   -2  1   0   -3  0   5   -2  -1  -3  -3  -4  
S   0   3   0   -1  -7  0   -2  -4  0   -5  -2  1   1   0   -2  4   3   -2  -5  -5  
T   0   0   0   -1  -6  -1  -1  -3  1   -5  -1  1   0   1   -1  3   3   -1  -5  -5  
V   -1  -2  -3  -3  -2  -6  -3  3   -3  -2  -2  -2  -1  -2  -3  -2  -1  6   -5  -5  
W   -5  -3  -4  -6  0   -5  -1  -7  -4  -4  -2  -3  -3  -3  -3  -5  -5  -5  10  -4  
Y   -7  2   -6  -7  3   -5  3   -5  -4  -4  -4  -4  -5  -4  -4  -5  -5  -5  -4  6   
</code>

Et en utilisant le BLOSUM 62:

<code python>
Substitution matrix generated with blocks from the Tyrosine kinase family
Identity% = 62.0%
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y   
A   2   1   0   1   -4  -1  3   -2  1   -2  -3  1   0   1   -1  0   0   -1  -14 0   
C   1   7   -10 -8  -12 -9  -10 -3  -8  -4  -7  -8  -3  -8  -10 -7  0   -3  -10 -12 
D   0   -10 4   2   -8  1   0   -5  2   -3  -3  3   -2  2   0   1   -2  -4  -24 -2  
E   1   -8  2   3   -5  2   1   -4  2   -2  -3  2   -3  2   0   2   -1  -3  -15 -1  
F   -4  -12 -8  -5  8   -7  -4  1   -6  -1  -2  -7  -12 -6  -7  -7  -7  0   -8  1   
G   -1  -9  1   2   -7  8   0   -6  1   -2  -3  1   -6  1   -1  2   2   -4  -9  -4  
H   3   -10 0   1   -4  0   4   -5  2   -1  -5  2   -4  2   0   0   0   -5  -8  2   
I   -2  -3  -5  -4  1   -6  -5  5   -4  4   1   -4  -2  -4  -6  -6  -4  5   -6  -2  
K   1   -8  2   2   -6  1   2   -4  2   0   -3  2   -5  2   1   1   1   -2  -11 -2  
L   -2  -4  -3  -2  -1  -2  -1  4   0   4   1   -1  -4  0   -2  -2  -1  4   -6  -1  
M   -3  -7  -3  -3  -2  -3  -5  1   -3  1   7   -3  -14 -2  -4  -3  -3  1   -8  -4  
N   1   -8  3   2   -7  1   2   -4  2   -1  -3  2   -5  2   0   1   0   -3  -14 -2  
P   0   -3  -2  -3  -12 -6  -4  -2  -5  -4  -14 -5  5   -4  -6  1   -1  -3  -20 -7  
Q   1   -8  2   2   -6  1   2   -4  2   0   -2  2   -4  2   1   1   -2  -3  -11 -1  
R   -1  -10 0   0   -7  -1  0   -6  1   -2  -4  0   -6  1   5   -1  -2  -5  -13 -2  
S   0   -7  1   2   -7  2   0   -6  1   -2  -3  1   1   1   -1  2   3   -2  -19 -3  
T   0   0   -2  -1  -7  2   0   -4  1   -1  -3  0   -1  -2  -2  3   6   0   -11 -4  
V   -1  -3  -4  -3  0   -4  -5  5   -2  4   1   -3  -3  -3  -5  -2  0   5   -4  -1  
W   -14 -10 -24 -15 -8  -9  -8  -6  -11 -6  -8  -14 -20 -11 -13 -19 -11 -4  9   -9  
Y   0   -12 -2  -1  1   -4  2   -2  -2  -1  -4  -2  -7  -1  -2  -3  -4  -1  -9  7   
Substitution matrix generated with blocks from the SH2 family
Identity% = 62.0%
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y   
A   5   -1  1   0   -5  0   -2  -3  0   -5  2   1   1   0   -1  0   0   -1  -4  -7  
C   -1  3   1   -3  1   -2  -2  0   0   -2  0   0   0   -1  0   2   0   -2  -1  1   
D   1   1   3   3   -7  1   -1  -3  1   -5  -1  2   1   2   -1  0   0   -3  -2  -6  
E   0   -3  3   5   -10 -1  -3  -2  2   -5  -1  0   -1  2   -1  -1  0   -4  -5  -8  
F   -5  1   -7  -10 7   -5  -1  -6  -4  -3  -2  -4  -6  -6  -4  -7  -6  -2  1   2   
G   0   -2  1   -1  -5  6   -5  -4  0   -5  -4  2   1   0   0   0   0   -6  -3  -5  
H   -2  -2  -1  -3  -1  -5  5   -6  -1  -5  0   0   -3  -1  -2  -2  -1  -3  -1  3   
I   -3  0   -3  -2  -6  -4  -6  5   -1  1   -3  -4  0   -2  -2  -4  -3  3   -6  -5  
K   0   0   1   2   -4  0   -1  -1  3   -4  0   1   0   2   2   0   1   -3  -3  -5  
L   -5  -2  -5  -5  -3  -5  -5  1   -4  5   3   -5  -1  -3  -2  -5  -5  -1  -3  -4  
M   2   0   -1  -1  -2  -4  0   -3  0   3   1   -1  0   0   1   -1  0   -2  -2  -3  
N   1   0   2   0   -4  2   0   -4  1   -5  -1  2   2   2   0   1   1   -2  -3  -4  
P   1   0   1   -1  -6  1   -3  0   0   -1  0   2   5   2   -2  1   0   -1  -4  -5  
Q   0   -1  2   2   -6  0   -1  -2  2   -3  0   2   2   2   1   0   1   -2  -2  -5  
R   -1  0   -1  -1  -4  0   -2  -2  2   -2  1   0   -2  1   5   -2  0   -3  -1  -3  
S   0   2   0   -1  -7  0   -2  -4  0   -5  -1  1   1   0   -2  4   3   -2  -5  -4  
T   0   0   0   0   -6  0   -1  -3  1   -5  0   1   0   1   0   3   3   -2  -4  -6  
V   -1  -2  -3  -4  -2  -6  -3  3   -3  -1  -2  -2  -1  -2  -3  -2  -2  6   -5  -5  
W   -4  -1  -2  -5  1   -3  -1  -6  -3  -3  -2  -3  -4  -2  -1  -5  -4  -5  10  -3  
Y   -7  1   -6  -8  2   -5  3   -5  -5  -4  -3  -4  -5  -5  -3  -4  -6  -5  -3  6
</code>
==== 3. Quelques Exemples d'alignement pour des séquences d'une même famille() ====
Y a t-il des différences entre les alignements quand vous utilisez les matrices 40% ou 70% ?

En utilisant l'algorithme de projet n°1, essayons de voir les différences entre l'utilisation des différentes matrices:

La famille utilisée est SH2

-Matrice 40%:

<code>
Welcome to the sequence aligner,
 Introduce the following information to obtain the alignment between 2 sequences.
Enter the 1st sequence filename: data MP1 2015/SH2-sequences.fasta
Choose 1 sequence from all given in the file ( from 0 to 5 ): 0
Same thing for the 2nd sequence,
Enter the 2nd sequence filename: data MP1 2015/SH2-sequences.fasta
Choose 1 sequence from all given in the file ( from 0 to 5 ): 1
Enter the filename of the score matrix you want to use: ScoreMatrix/blosum40.txt
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
score : 257
</code>

-Matrice 70%: 

<code>
Welcome to the sequence aligner,
 Introduce the following information to obtain the alignment between 2 sequences.
Enter the 1st sequence filename: data MP1 2015/SH2-sequences.fasta
Choose 1 sequence from all given in the file ( from 0 to 5 ): 0
Same thing for the 2nd sequence,
Enter the 2nd sequence filename: data MP1 2015/SH2-sequences.fasta
Choose 1 sequence from all given in the file ( from 0 to 5 ): 1
Enter the filename of the score matrix you want to use: ScoreMatrix/blosum70.txt            
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
score : 311
</code>

==== 4. Comparaison des alignements pour les mêmes séquences en utilisant BLOSUM62 ====

Les alignements obtenus en utilisant les matrices construites sont-ils meilleurs ?

On va utiliser différentes matrices créés et les comparer aux matrices données pour voir la différence:

- Substitution matrix generated with blocks from the Tyrosine kinase family BLOUM62:

<code>
Welcome to the sequence aligner,
 Introduce the following information to obtain the alignment between 2 sequences.
Enter the 1st sequence filename: data MP1 2015/SH2-sequences.fasta
Choose 1 sequence from all given in the file ( from 0 to 5 ): 0
Same thing for the 2nd sequence,
Enter the 2nd sequence filename: data MP1 2015/SH2-sequences.fasta
Choose 1 sequence from all given in the file ( from 0 to 5 ): 1
Enter the filename of the score matrix you want to use: ScoreMatrix/blosum62TK.txt
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
score : 316
</code>


- Substitution matrix generated with blocks from the SH2 family BLOSUM62:

<code>
Welcome to the sequence aligner,
 Introduce the following information to obtain the alignment between 2 sequences.
Enter the 1st sequence filename: data MP1 2015/SH2-sequences.fasta
Choose 1 sequence from all given in the file ( from 0 to 5 ): 0
Same thing for the 2nd sequence,
Enter the 2nd sequence filename: data MP1 2015/SH2-sequences.fasta
Choose 1 sequence from all given in the file ( from 0 to 5 ): 1
Enter the filename of the score matrix you want to use: ScoreMatrix/blosum62SH.txt
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
score : 315
</code>

- Given Matrix BLOSUM62:
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

On peut clairement remarquer que nos matrices sont très similaires dans leur effet mais que celle qui nous a été fournie est bien meilleure.