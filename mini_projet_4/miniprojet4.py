from math import sqrt,log

def parsageDSSP(nomFichier):
	liste = []
	try:
		fichier = open("dataset/"+nomFichier,"r")
		for line in fichier:
			liste.append(line.strip().split())
		for i in range(len(liste)):
			liste[i] = [liste[i][0][0:4],liste[i][0][4]]
		#return une liste de string avec le nom du dssp à utiliser ainsi que la chaîne à extraire
		return liste
	except(IOError):
		print("Erreur lors de l'ouverture du fichier '"+nomFichier+"'")
		return None

def traitementDSSP(listeDSSP,nomFichier):
	try:
		if(nomFichier == "CATH_info_test.txt"):
			fichierSequence = "sequencesTest.txt"
			dossier = "dataset/dssp_test/"
		elif(nomFichier == "CATH_info.txt"):
			fichierSequence = "sequencesEntrainement.txt"
			dossier = "dataset/dssp/"
		#liste[i] = [nomDSSP,chaine]
		fichier = open(fichierSequence,"w")
		for i in range(len(listeDSSP)):
			nomProtein,organisme,sequence,structureSecondaire = lectureDonnees(dossier,listeDSSP[i][0],listeDSSP[i][1])
			fichier.write("> "+listeDSSP[i][0]+" | "+nomProtein+" | "+organisme+" \n")
			fichier.write(sequence.upper())
			fichier.write("\n")
			fichier.write(structureSecondaire)
			fichier.write("\n")
		fichier.close()		
	except(IOError):
		print("Erreur lors du traitement des dssp")

def lectureDonnees(dossier,fichierDSSP,chaine):
	try:
		fichier = open(dossier+fichierDSSP+".dssp","r")
		lignes = fichier.readlines()
		organisme=""
		sequence = ""
		secondaire = ""
		nomProtein = lignes[3].split()[3]
		if(nomProtein[-1] == ";" or nomProtein[-1] == ","):
			nomProtein = nomProtein[:-1]
		listeOrganisme = lignes[4].split()[3:]
		for mot in listeOrganisme:
			if(mot[-1] == "." or mot[-1] == ";"):
				mot = mot[:-1]
			organisme+=(mot+" ")
		for i in range(28,len(lignes)):
			if(lignes[i][11] == chaine):
				sequence+=lignes[i][13]
				if(lignes[i][16] == "G" or lignes[i][16] == "I"):
					secondaire+="H"
				elif(lignes[i][16] == "S" or lignes[i][16] == "B" or lignes[i][16] == " "):
					secondaire+="C"
				else:
					secondaire+=lignes[i][16]
		fichier.close()
		return nomProtein,organisme,sequence.upper(),secondaire
	except(IOError):
		print("Erreur lors de l'ouverture de "+nomFichier+".dssp")
		return None

def creationDico(dico,nomFichier):
	listeSequences,listeStructures = lectureListes(nomFichier)
	for i in range(len(listeSequences)):
		for j in range(len(listeSequences[i])):
			if ((listeStructures[i][j]) not in dico):
				dico[(listeStructures[i][j])] = 1
			else:
				dico[(listeStructures[i][j])] += 1		
			for k in range (max(0, j - 8), min(len(listeSequences[i]), j + 8)):
				if (j == k):
					if ((listeStructures[i][j], listeSequences[i][k]) not in dico):
						dico[(listeStructures[i][j], listeSequences[i][j])] = 1
					else:
						dico[(listeStructures[i][j], listeSequences[i][j])] += 1
				else:
					if ((listeStructures[i][j], listeSequences[i][k], listeSequences[i][j]) not in dico):
						dico[(listeStructures[i][j], listeSequences[i][k], listeSequences[i][j])] = 1
					else:
						dico[(listeStructures[i][j], listeSequences[i][k], listeSequences[i][j])] += 1


def lectureListes(nomFichier):
	listeSequences = [] 
	listeStructures = []
	try:
		fichier = open(nomFichier,"r")
		donnees = []
		lecture = ""
		for ligne in fichier:
			if(ligne[0] == ">"):
				donnees.append(lecture.split())
				lecture = ""
			else:
				lecture+=ligne
		donnees.append(lecture.split())
		for i in range(1,len(donnees)):
			listeSequences.append(donnees[i][0])
			listeStructures.append(donnees[i][1])
		fichier.close()
		return listeSequences, listeStructures
	except (IOError):
		# Le fichier n'existe pas ou permission non accordée
		print ("Erreur lors de l'ouverture du fichier!")
		return None

def notStructure(dico, structure):
	occurenceH = dico['H']
	occurenceE = dico['E']
	occurenceC = dico['C']
	occurenceT = dico['T']  
	if (structure == 'H'):
		resultat = occurenceE + occurenceC + occurenceT
	elif (structure == 'E'):
		resultat = occurenceH + occurenceC + occurenceT
	elif (structure == 'C'):
		resultat = occurenceE + occurenceH + occurenceT
	elif (structure == 'T'):
		resultat = occurenceE + occurenceC + occurenceH   
	return resultat

def notStructureAcid(dico, structure, acid):
	occurenceH_Acid = dico[('H', acid.upper())]
	occurenceE_Acid = dico[('E', acid.upper())]
	occurenceC_Acid = dico[('C', acid.upper())]
	occurenceT_Acid = dico[('T', acid.upper())]
	if (structure == 'H'):
		resultat = occurenceE_Acid + occurenceC_Acid + occurenceT_Acid
	elif (structure == 'E'):
		resultat = occurenceH_Acid + occurenceC_Acid + occurenceT_Acid
	elif (structure == 'C'):
		resultat = occurenceE_Acid + occurenceH_Acid + occurenceT_Acid
	elif (structure == 'T'):
		resultat = occurenceE_Acid + occurenceC_Acid + occurenceH_Acid
	return resultat
	
def notStructureAcid1Acid2(dico, structure, acid2, acid1):  
	occurenceH_Acid1_Acid2 = dico[('H', acid2.upper(), acid1.upper())]
	occurenceE_Acid1_Acid2 = dico[('E', acid2.upper(), acid1.upper())]
	occurenceC_Acid1_Acid2 = dico[('C', acid2.upper(), acid1.upper())]
	occurenceT_Acid1_Acid2 = dico[('T', acid2.upper(), acid1.upper())]
	if (structure == 'H'):
		resultat = occurenceE_Acid1_Acid2 + occurenceC_Acid1_Acid2 + occurenceT_Acid1_Acid2
	elif (structure == 'E'):
		resultat = occurenceH_Acid1_Acid2 + occurenceC_Acid1_Acid2 + occurenceT_Acid1_Acid2
	elif (structure == 'C'):
		resultat = occurenceE_Acid1_Acid2 + occurenceH_Acid1_Acid2 + occurenceT_Acid1_Acid2
	elif (structure == 'T'):
		resultat = occurenceE_Acid1_Acid2 + occurenceC_Acid1_Acid2 + occurenceH_Acid1_Acid2
	return resultat
	
def lePlusGrand(dico, index):
	liste = []
	listeStructures = []
	for key in dico:
		if (key[1] == index):
			liste.append(dico[key])
			listeStructures.append(key[0])
	maximum = -10
	for index in range (len(liste)):
		if (liste[index] > maximum):
			maximum = liste[index]
			i = index
	return listeStructures[i]

def informationsMutuelles(dico, dicoInfo, sequence, structure):
	#parcours sequence
	for i in range (len(sequence)):
		infoIndividuelle = 0
		AA1 = sequence[i]
		terme1 = log(dico[(structure, AA1.upper())] / (notStructureAcid(dico, structure, AA1.upper())))
		terme2 = log((notStructure(dico, structure)) / dico[structure])
		infoIndividuelle = terme1 + terme2
		infoPaires = 0
		for j in range (max(0, i - 8), min(len(sequence), i + 8)):
			if (j != 0):
				AA2 = sequence[j]
				terme3 = log(dico[(structure, AA2.upper(), AA1.upper())] / (notStructureAcid1Acid2(dico, structure, AA2.upper(), AA1.upper())))
				terme4 = log((notStructureAcid(dico, structure, AA1.upper())) / dico[(structure, AA1.upper())])
				infoPaires += (terme3 + terme4)
		dicoInfo[(structure, i)] = infoIndividuelle + infoPaires

def resultatComparaison(bonneStructure, prediction, forme):
	truePositive = 0
	trueNegative = 0
	falsePositive = 0
	falseNegative = 0
	for i in range (len(bonneStructure)):
		if (bonneStructure[i] == prediction[i] == forme):
			truePositive += 1
		elif (bonneStructure[i] == prediction[i] != forme):
			trueNegative += 1
		elif (bonneStructure[i] != prediction[i]) and (prediction[i] == forme):
			falsePositive += 1
		elif (bonneStructure[i] == forme) and (prediction[i] != forme):
			falseNegative += 1	
	return truePositive, trueNegative, falsePositive, falseNegative

def comparaisonListes(bonneStructure, prediction):
	compteur = 0
	for i in range (len(bonneStructure)):
		if (bonneStructure[i] == prediction[i]):
			compteur += 1 
	return compteur

def mathewsCorrelationCoefficient(truePositive, trueNegative, falsePositive, falseNegative):
	try:
		num = (truePositive * trueNegative) - (falsePositive * falseNegative)
		denum = sqrt((truePositive + falsePositive) * (truePositive + falseNegative) 
				  * (trueNegative + falsePositive) * (trueNegative + falseNegative))
		MCC = num/denum 
		return MCC        
	except (ZeroDivisionError):
		return None

def Q3(predictionCorrecte,nombreResidus):
	return predictionCorrecte/nombreResidus*100

def GORIII(dico,nomFichier):
	listeSequences,listeStructures = lectureListes(nomFichier)
	listeFormes = ['C', 'E', 'H', 'T']
	for i in range(len(listeSequences)):
		infosMutuelles = {}
		prediction = ""
		for forme in listeFormes:
			informationsMutuelles(dico,infosMutuelles,listeSequences[i],forme)
		for j in range(len(listeSequences[i])):
			prediction+= lePlusGrand(infosMutuelles,j)
		predictionCorrecte = comparaisonListes(listeStructures[i],prediction)
		q3 = Q3(predictionCorrecte,len(listeSequences[i]))
		valeurs = []
		for forme in listeFormes:
			truePositive, trueNegative, falsePositive, falseNegative = resultatComparaison(listeStructures[i], prediction, forme)
			MCC = mathewsCorrelationCoefficient(truePositive, trueNegative, falsePositive, falseNegative)
			valeurs.append(MCC)
		print ("Q3_GLOBAL =", q3, "MCC_C =", valeurs[0], "MCC_E =", valeurs[1], "MCC_H =", valeurs[2], "MCC_T =", valeurs[3])
				


if (__name__ == "__main__"):
	print("Bonjour et bienvenu dans ce magnifique programme, j'espère qu'il fait beau et que vous profitez du soleil mais pas trop. Pour en revenir à la partie intéressante, vous avez plus ou moins 15 secondes pour faire autre chose pendant que le programme s'exécute (ex: lire le journal, penser à l'avenir du monde, à limmensité de l'univers, bref vous avez le choix")
	dico = {}
	print("Création d'une liste d'objets DSSP grâce au fichier 'CATH_info.txt'...")
	listeDSSP = parsageDSSP("CATH_info.txt")
	print("Traitement des objets DSSP avec les données d'entraînement...")
	traitementDSSP(listeDSSP,"CATH_info.txt")
	print("Création d'une liste d'objets DSSP grâce au fichier 'CATH_info_test.txt'...")
	listeDSSP = parsageDSSP("CATH_info_test.txt")
	print("Traitement des objets DSSP du fichier test...")
	traitementDSSP(listeDSSP,"CATH_info_test.txt")
	print("Création du dictionnaire nécéssaire au calcul des fréquences et utilisé par l'algorithme GOR III ...")
	creationDico(dico,"sequencesEntrainement.txt")
	print("Application de l'algorithme GOR III afin de prédire la structure secondaire des séquences de test...")
	GORIII(dico,"sequencesTest.txt")