from Alignment import *


def main():
	print("Welcome to the sequence aligner,\n Introduce the following information to obtain the alignment between 2 sequences.")

	sequenceFilename = str(input("Enter the 1st sequence filename: "))
	sequences1 = Sequence.extraction(sequenceFilename)
	sequenceIndex = int(input("Choose 1 sequence from all given in the file ( from 0 to "+str(len(sequences1)-1)+" ): "))
	sequence1 = sequences1[sequenceIndex]


	sequenceFilename = str(input("Same thing for the 2nd sequence,\nEnter the 2nd sequence filename: "))
	sequences2 = Sequence.extraction(sequenceFilename)
	sequenceIndex = int(input("Choose 1 sequence from all given in the file ( from 0 to "+str(len(sequences2)-1)+" ): "))
	sequence2 = sequences2[sequenceIndex]


	scoreFilename = str(input("Enter the filename of the score matrix you want to use: "))
	score = Score.extraction(scoreFilename)
	

	gapStart = -1
	while gapStart < 0:
		gapStart = int(input("Enter the opening gap penalty: "))

	gapExtend = -1
	while gapExtend < 0:
		gapExtend = int(input("Enter the extending gap penalty: "))


	alignmentType = ""
	while not(alignmentType == "local" or alignmentType == "global"):
		alignmentType = str(input("Give the alignmentType, 'global' or 'local': "))


	alignment = Alignment(sequence1,sequence2,score,gapStart,gapExtend,alignmentType)
	alignment.align()
	print(alignment)

if __name__ == '__main__':
	main()
	