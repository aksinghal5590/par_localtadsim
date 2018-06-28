import sys

def countHangingTADs(inFile1, inFile2, outFile):

	factor = 100000

	inList1 = list()
	inList2 = list()
	outList = list()

	for line in open(inFile1):

		data = line.split()
		inList1.append((int(data[1])/factor, (int(data[2]) + 1)/factor))

	for line in open(inFile2):
		data = line.split()
		inList2.append((int(data[1])/factor, (int(data[2]) + 1)/factor))

	p = 0
	for line in open(outFile):
		if p < 1:
			p = p + 1
			continue
		data = line.split()
		outList.append((int(data[0]), int(data[1])))

	htads = list()
	for tad in outList:
		for in_tad in inList1:
			if in_tad[0] < tad[0] < in_tad[1]:
				if tad[0] > (in_tad[1] - in_tad[0])/2 + in_tad[0]:
					htads.append(tad)
			if in_tad[0] < tad[1] < in_tad[1]:
				if tad[1] < (in_tad[1] - in_tad[0])/2 + in_tad[0]:
					htads.append(tad)

		for in_tad in inList2:
			if in_tad[0] < tad[0] < in_tad[1]:
				if tad[0] > (in_tad[1] - in_tad[0])/2 + in_tad[0]:
					htads.append(tad)
			if in_tad[0] < tad[1] < in_tad[1]:
				if tad[1] < (in_tad[1] - in_tad[0])/2 + in_tad[0]:
					htads.append(tad)

	print(list(set(htads)))



countHangingTADs(sys.argv[1], sys.argv[2], sys.argv[3])