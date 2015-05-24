#On place toutes les séquences téléchargez dans des blocks différents
def divisionBlocks(nomFichier):
	blockA = []
	blockB = []
	blockC = []
	blockD = []
	blockE = []

	listeBlocks = [blockA,blockB,blockC,blockD,blockE]

	try:
		fichier = open(nomFichier,"r")
		compteur = 0
		sequence = " "
		while(sequence != ""):
			sequence = fichier.readline().strip("\n")
			if(compteur%6 != 0):
				listeBlocks[compteur%6-1].append(sequence)
			compteur+=1
		fichier.close()
		return listeBlocks
	except(IOError):
		print("IOError")
		return None

#On crée un fichir pour chaque block différent
def ecritureBlocks(listeBlocks, nomFichierDeBase):
	print("NOM: "+ nomFichierDeBase)
	indice = 0
	listeSuffixes = ["A","B","C","D","E"]
	for block in listeBlocks:
		try:
			nomFichierCree = nomFichierDeBase+listeSuffixes[indice]+".txt"
			indice+=1 
			fichier = open(nomFichierCree,"w")
			for sequence in block:
				fichier.write(sequence+"\n")
			fichier.close()
		except(IOError):
			print("IOError")

if(__name__ == "__main__"):
	ecritureBlocks(divisionBlocks("PR00109.txt"), "PR00109")
	ecritureBlocks(divisionBlocks("PR00401.txt"), "PR00401")
