from math import sqrt,log

def parser(filename):
	listeSequences = []
	try:
		fichier = open(filename,"r")
		sequence = ""
		for ligne in fichier:
			if(ligne[0] != ">"):
				sequence+=ligne.strip("\n")
			else:
				if(sequence != ""):
					listeSequences.append(sequence)
					sequence = ""
		return listeSequences
	except IOError as e:
		print(e)
		return None

def calculAA(liste):
	listeAcidesAmines = []
	for sequence in liste:
		for aa in sequence:
			if(aa not in listeAcidesAmines and aa != "-"):
				listeAcidesAmines.append(aa)
	return listeAcidesAmines

def MatriceAAPosition(liste):
	listeAA = calculAA(liste)
	mat = [[0 for i in range(len(liste[0]))] for j in range(len(listeAA))]
	for sequence in (liste):
		for indicePosition,aa in enumerate(sequence):
			if( aa != "-"):
				mat[listeAA.index(aa)][indicePosition]+=1
	return mat

def calculDesFrequences(matrice,nombreSequences):
	for i in range(len(matrice)):
		for j in range(len(matrice[0])):
			matrice[i][j] /= nombreSequences
	return matrice

def ajoutPseudoCounts(matrice,nombreSequences,listeSwissProt,listeAA):
	for i in range(len(matrice)):
		for j in range(len(matrice[0])):
			matrice[i][j] = (nombreSequences - 1)*(matrice[i][j]) + (sqrt(nombreSequences))*(listeSwissProt[listeAA[i]])/(sqrt(nombreSequences)+(nombreSequences - 1))
	return matrice		

def matrice_M(matrice,listeSwissProt,listeAA):
	for i in range(len(matrice)):
		for j in range(len(matrice[0])):
			if(listeSwissProt[listeAA[i]] == 0):
				matrice[i][j] = 0
			else:
				matrice[i][j] = log(matrice[i][j]/(listeSwissProt[listeAA[i]]/100),10)
	return matrice

def consensus(matrice,listeAA):
	consensus = []
	print(matrice[0][0])
	for i in range(len(matrice[0])):
		maxi = [0,0] # [score,position]
		for j in range(len(matrice)):
			if(matrice[j][i] > maxi[0]):
				maxi[0] = matrice[j][i]
				maxi[1] = j
		consensus.append(listeAA[maxi[1]])
	return consensus


def writeToFile(matrix, fileName,listeAA):
    try:
        file = open(fileName, 'w')
        lineIndex = 0
        while (lineIndex != len(matrix)):
            columnIndex = 0
            while (columnIndex != len(matrix[0])):
                if (columnIndex == 0) and (lineIndex == 0):
                    file.write('- ')
                elif (columnIndex == 0):
                    file.write('\n')
                    file.write(listeAA[lineIndex - 1])
                    file.write(' ')
                else:
                    file.write(str(float(matrix[lineIndex][columnIndex])))
                    file.write(' ')
                columnIndex = columnIndex + 1
            lineIndex = lineIndex + 1
            file.write(' ')
        file.close()
        print ("Ecriture terminée")
    except (IOError):
        print ("Impossible d'écrire dans le fichier")
        return None

if(__name__ == "__main__"):
	listeSwissProt = {"A": 8.28, "Q": 3.94, "L": 9.67, "S": 6.50, "R": 5.53, "E": 6.76, "K": 5.85, "T": 5.32, "N": 4.05, "G": 7.09, "M": 2.43, "W": 1.07, "D": 5.45, "H": 2.27, "F": 3.86, "Y": 2.91, "C": 1.36, "I": 5.99, "P": 4.68, "V": 6.87, "B": 0, "Z": 0, "X": 0}
	quit = False
	while( not quit):
		choix = input("Bonjour, bienvenu dans le créateur de profils pour des séquences de même famille. Veuillez rentrer l'outil avec lequel vous désirez aligner les liste de séquences:\nSH2\n   (1.1)MUSCLE\n   (1.2)CLUSTAL\nSH3\n   (2.1)MUSCLE\n   (2.2)CLUSTAL\nChoix(ex '1.1' pour le domaine SH2 aligné avec l'outil MUSCLE): ")
		if choix == '1.1':
			listeDeSequences = parser("msaresults-MUSCLE.fasta")
			titre = "PSSM_SH2_MUSCLE.txt"
		elif choix == '1.2':
			listeDeSequences = parser("msaresults-CLUSTAL.fasta")
			titre = "PSSM_SH2_CLUSTAL.txt"
		elif choix == '2.1':
			listeDeSequences = parser("msaresults2-MUSCLE.fasta")
			titre = "PSSM_SH3_MUSCLE.txt"
		elif choix == '2.2':
			listeDeSequences = parser("msaresults2-CLUSTAL.fasta")
			titre = "PSSM_SH3_CLUSTAL.txt"
		else:
			quit = True
			print("A bientôt sur ce magnifique projet.")
		if(quit == False):
			matriceDeBase = MatriceAAPosition(listeDeSequences)
			matriceDesFrequences = calculDesFrequences(matriceDeBase,len(listeDeSequences))
			matricePseudoCount = ajoutPseudoCounts(matriceDesFrequences,len(listeDeSequences),listeSwissProt,calculAA(listeDeSequences))
			PSSM = matrice_M(matricePseudoCount,listeSwissProt,calculAA(listeDeSequences))
			verification = consensus(PSSM,calculAA(listeDeSequences))
			print(len(PSSM[0] ))
			writeToFile(PSSM,titre,calculAA(listeDeSequences))
			bornes = input("Veuillez rentrer les bornes des positions que vous désirez comparer avec le weblogo (ex: '5 10'): ")
			start,end = bornes.split()
			print(int(start),int(end),len(verification))
			print(int(start) < 0)
			print(int(end) > len(verification))
			while(int(start) < 0 or int(end) > len(verification)):
				print("Vous avez dépassé les bornes du consensus(1:"+str(len(verification))+")\n")
				bornes = input("Veuillez rentrer les bornes des positions que vous désirez comparer avec le weblogo (ex: '5 10'): ")
				start,end = bornes.split()
			print(verification[(int(start)-1):int(end)]) #indice weblogo commence à 1 et la list à 0 => weblogo de 890 à 900
