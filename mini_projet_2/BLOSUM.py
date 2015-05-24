import copy
from math import log

def BLOSUM(listeFichiers,pourcentage):
	listeMatrices = []
	for nomFichier in listeFichiers:
		#######listeSequences#######
		listeSequences = recupererSequences(nomFichier)
		#######groupesDifferents#######
		listeGroupes = trouverGroupesDifferents(listeSequences,pourcentage)
		#######ajoutListeMatrices#######
		listeMatrices.append(matriceFrequencesPonderees(listeGroupes,['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']))
	#######moyenneDesMatrices#######
	matriceResultante = moyenneMatrices(listeMatrices)
	#######probabilitésD'Occurences#######
	matriceOccurences = calculProbabiliteOccurence(matriceResultante)
	#afficherMatrice(matriceOccurences)
	#######calculDeLogChance#######
	matriceFinale = calculMatriceFinale(matriceOccurences)
	afficherMatrice(matriceFinale)
	

def calculMatriceFinale(matriceOccurences):
	copieMatrice = copy.deepcopy(matriceOccurences)
	for i in range(len(matriceOccurences)):
		for j in range(len(matriceOccurences[0])):
			matriceOccurences[i][j] = calculTauxLogChance(copieMatrice,i,j)
	return matriceOccurences


def calculTauxLogChance(matriceOccurences,i,j):
	try:
		#print("FREQUENCE")
		#print(frequencePrevuePourAlignement(matriceOccurences,i,j))
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
	
	
def indiceProteine(proteine):
	listeProteines = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	indice = 0
	for i in range(len(listeProteines)):
		if(proteine == listeProteines[i]):
			indice = i
	return indice

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


#Moyenne des matrices passées en paramètres
def moyenneMatrices(listeMatrices):
	matrice = [[0 for i in range(20)] for j in range(20)]
	for i in range(len(listeMatrices[0])):
		for j in range(len(listeMatrices[0][0])):
			for n in range(len(listeMatrices)):

				matrice[i][j] += listeMatrices[n][i][j]
			matrice[i][j] = matrice[i][j]/len(listeMatrices)
	return matrice

#Récupère les séquences dans un fichier donné
def recupererSequences(nomFichier):
	try:
		fichier = open(nomFichier,"r")
		listeSequences = fichier.readlines()
		for i in range(len(listeSequences)):
			listeSequences[i] = listeSequences[i].strip("\n")
		fichier.close()
		return listeSequences
	except (IOError):
		print ("IOError") 
		return None

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



def listeProteines(listeSequences):
    liste = []
    for sequence in listeSequences:
        for proteine in sequence:
            if (proteine not in liste):
                liste.append(proteine)
    return liste


def afficherMatrice(matrice):
	listeProteines = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	print(end='   ')
	for proteine in listeProteines:
		print(proteine,end = '   ')
	print()
	for i in range(len(matrice)):
		print(listeProteines[i], end = '   ')
		for j in range(len(matrice[0])):
			print(round(matrice[i][j]),end = ' '*(4-len(str(round(matrice[i][j])))))
		print()


if(__name__ == "__main__"):
	print("Substitution matrix generated with blocks from the Tyrosine kinase family")
	print("Identity% = 40.0%")
	BLOSUM(["PR00109A.txt","PR00109B.txt","PR00109C.txt","PR00109D.txt","PR00109E.txt"],0.4)
	print("Substitution matrix generated with blocks from the Tyrosine kinase family")
	print("Identity% = 70.0%")	
	BLOSUM(["PR00109A.txt","PR00109B.txt","PR00109C.txt","PR00109D.txt","PR00109E.txt"],0.7)
	print("Substitution matrix generated with blocks from the SH2 family")
	print("Identity% = 40.0%")	
	BLOSUM(["PR00401A.txt","PR00401B.txt","PR00401C.txt","PR00401D.txt","PR00401E.txt"],0.4)
	print("Substitution matrix generated with blocks from the SH2 family")
	print("Identity% = 70.0%")	
	BLOSUM(["PR00401A.txt","PR00401B.txt","PR00401C.txt","PR00401D.txt","PR00401E.txt"],0.7)
