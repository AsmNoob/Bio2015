from Score import *
from Sequence import *


class Alignment:
	def __init__(self, sequence1,sequence2, score, gapStart, gapExtend, alignmentType):
		self.sequence1 = sequence1
		self.sequence2 = sequence2
		self.score = score
		self.gapStart = gapStart
		self.gapExtend = gapExtend
		self.alignmentType = alignmentType
		self.lengths = [len(sequence1),len(sequence2)]
		self.finalAlignment = []


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


	def isFromUp(self,i,j):
		return (j > 0 and self.S[i][j] == self.W[i][j])


	def isFromDiagonal(self,i,j):
		return (i > 0 and j > 0 and self.S[i][j] == self.S[i-1][j-1]+self.score[self.sequence1[i-1],self.sequence2[j-1]])


	def isFromLeft(self,i,j):
		return (i > 0 and self.S[i][j] == self.V[i][j])


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

	def align(self):
		self.fillMatrix()
		self.traceback(self.lengths[0],self.lengths[1],"","")

	def __repr__(self):
		res = ""

		for i in range(len(self.finalAlignment)):
			res += self.finalAlignment[i][0] + "\n"
			identity = 0
			similarity = 0
			gap = 0
			for j in range(len(self.finalAlignment[i][0])):
				if (self.finalAlignment[i][0][j] == self.finalAlignment[i][1][j]):
					res += ":"
					identity += 1
					similarity += 1	
				elif (self.finalAlignment[i][0][j] == '-' or self.finalAlignment[i][1][j] == '-'):
					res += " "
					gap += 1	
				else:
					res += "."
					similarity += 1
			res += "\n" + self.finalAlignment[i][1] + "\n"
			res += str(round(100*identity/len(self.finalAlignment[i][0]), 1)) + "% identity\n"
			res += str(round(100*similarity/len(self.finalAlignment[i][0]), 1)) + "% similarity\n"
			res += str(round(100*gap/len(self.finalAlignment[i][0]), 1)) + "% gap\n"
			res += "Length : " + str(len(self.finalAlignment[i][0])) + "\n"
			res += "score : " + str(self.S[self.lengths[0]][self.lengths[1]])
			return res