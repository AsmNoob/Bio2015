class Score:
	def __init__(self, matrix, indexes):
		self.matrix = matrix
		self.indexes = indexes


	def extraction(filename):
		matrix = []
		indexes = ""
		try:
			file = open(filename,"r")
			
			for line in file:
				if line[0] != "#":
					if line[0] == " " or line[0] == "\t":
						indexes = line.split()
					elif line[0].isalpha() or line[0] == "*" or line[0] == "-":
						matrix.append(list(map(int, line[1:].split())))
		except FileNotFoundError as e:
			print("FileNotFoundError ({0}): {1}".format(e.errno,e.strerror))
		else:
			file.close()
		
		return Score(matrix,indexes)


	def __getitem__(self,acid):
		return self.matrix[self.indexes.index(acid[0])][self.indexes.index(acid[1])]


	def __repr__(self):
		res = " "
		for index in self.indexes:
			res+= "  " + index
		res+="\n"
		compteur = 0
		for line in self.matrix:
			res+= self.indexes[compteur]
			compteur+=1
			for char in line:
				if char >= 0:
					if char < 10:
						res+= " "
					res+= " " + str(char)
				elif char < 0:
					if abs(char) < 10 :
						res+= " "
					res+= str(char)
			res+= "\n"
		return res