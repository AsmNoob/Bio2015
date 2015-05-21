class Sequence(str):
	def __new__(cls, *args, **kw):
		return str.__new__(cls, *args, **kw)


	def extraction(filename):
		sequencesList = []
		sequence = ""
		try:
			file = open(filename,"r")
			for line in file:
				if line[0] != ">":
					sequence+= line[:-1].strip()
				else:
					if sequence != "":
						sequencesList.append(Sequence(sequence))
					sequence = ""
			if sequence != "":
				sequencesList.append(Sequence(sequence))
		except FileNotFoundError as e:
			print("FileNotFoundError ({0}): {1}".format(e.errno,e.strerror))
		else:
			file.close()
			return sequencesList


	def identify(self,other):
		res = 0
		if len(self) == len(other):
			for i in range(len(self)):
				if(self[i] == other[i]):
					res+=1
			res/= len(self)
		return res