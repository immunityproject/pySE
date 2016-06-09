from proteins import all_aminos,chains,chain_bases

from argparse import ArgumentParser


DEBUG = False

class Result:
	def __init__(self, proteinName, chain, offset, sequence, score):
		self.proteinName = proteinName
		self.chain = chain
		self.offset = offset
		self.sequence = sequence
		self.score = score
	

def getScore(seq1, seq2):
	
	total_match = 0.0

	if(len(seq1) != len(seq2)):
		print "Error, seq len mismatch:" + seq1 + ":" + seq2
	else:
		for i in range(len(seq1)):
			if (seq1[i-1] == seq2[i-1]):
				total_match = total_match + 1

	return total_match/len(seq1)

def scanSingleProteinForMatch(protein,proteinArray, epitope, cutoff):
	#print "Scanning protein " + protein
	resultArray = []

	#def __init__(self, proteinName, offset, sequence, score):
	#where we are, stop before len(epitope) so we don't go off the end of the amino array)

	for chain in proteinArray[protein]:
		chain_base = chain_bases[protein][chain] 
		for startLoc in range(len(proteinArray[protein][chain]) - len(epitope) + 1):
			testEpitope = proteinArray[protein][chain][startLoc:startLoc + len(epitope)]
			scoreResult = getScore(epitope,testEpitope)
			if (scoreResult >= cutoff):
				resultArray.append(Result(protein,chain,startLoc+chain_base,testEpitope,scoreResult))
			
			if(DEBUG):
				print "Added " + protein +"-"+chain + "-" + str(startLoc) + " (" + testEpitope + ") : " + str(scoreResult)	

	return resultArray
	

def scanProteinsForMatch(proteinArray, epitope, cutoff):
	
	resultsArray = []

	for protein in proteinArray:
		proteinResult = scanSingleProteinForMatch(protein, proteinArray, epitope, cutoff)
		resultsArray.extend(proteinResult)	

	resultsArray.sort(key = lambda result: 1-result.score)
	return resultsArray



def printTopResults(resultsArray, eptiope, numberOfResults):
	print "Best " + str(numberOfResults) + " matches for:" + eptiope
	
	for i in range(numberOfResults):
		if(i >= len(resultsArray)):
			print str(i) + ": ---"
		else:
			#def __init__(self, proteinName, offset, sequence, score):
			result = resultsArray[i]
			print str(i) + ": "  + result.proteinName + "-" + result.chain  + "-" + str(result.offset),
			print " (" + result.sequence +  ") " + "SCORE:" + str(result.score)  



def findChainsFromBestHit(bestResult):
	for chain in chains[bestResult.proteinName]:
		if bestResult.chain in chain:
			return chain


def constructResultDict(first,last,chain_array,conf,ref):
	return_dict = {}
	return_dict['first']=first
	return_dict['last']=last
	return_dict['chains']=chain_array
	return_dict['conf'] = conf
	return_dict['reference_value']=float(ref)
	
	return return_dict

def findLocalMatchFromList(epitope_array):

	#initialize protein dict.		
	protein_result = {};
	for protein in all_aminos:
		protein_result[protein] = {}

	
	for epitope_data in epitope_array:
		epitope = epitope_data[0]
		shortname = epitope[0] + epitope[-1] + str(len(epitope))
		resultsArray = scanProteinsForMatch(all_aminos, epitope, 0.2)
		#use best result and look up chain data.
		bestResult = resultsArray[0]
		chain_array = findChainsFromBestHit(bestResult)
		
		if bestResult.score < .665:
			print "Bad conf for: " + shortname + " " + str(bestResult.score)

		else:
			protein_result[bestResult.proteinName][shortname] = constructResultDict(bestResult.offset, bestResult.offset+len(bestResult.sequence)-1, chain_array, bestResult.score, epitope_data[1] )

	print protein_result

	#for protein in protein_result:
	#	for epitope in protein_result[protein]:
	#		print protein + "\t" + str(protein_result[protein][epitope]["chains"])
		

if __name__ == "__main__":

	parser = ArgumentParser()
	parser.add_argument('-e',dest='epitope',default=None, help = "list an epitope chain to search for. ex ALFALFQ")
	parser.add_argument('-c',dest='cutoff',default=None, help = "cutoff for if an match is interesting")
	parser.add_argument('-f',dest='filename',default=None, help = "specify a file to open that contains newline separated list of epitopes to search for")
	args = parser.parse_args()


	if args.epitope is not None:
		args.epitope = args.epitope.upper()
		
		print "I will search for: " + args.epitope

		print "I have the following proteins loaded into memory:"
		for protein in all_aminos:
			print "\t*" + protein
		
		resultsArray = scanProteinsForMatch(all_aminos, args.epitope, 0.2)
		printTopResults(resultsArray, args.epitope, 10)

	elif args.filename is not None:
		f=open(args.filename, 'r')
		epitope_array = []
		for line in f:
			split_line=line.strip().split(',')
			epitope_array.append([split_line[0],split_line[1]])

		findLocalMatchFromList(epitope_array)
	
			
		#except Exception,e:		
		#	print "Error opening " + args.filename
		#	print str(e)
	else:
		print "Use -e or -f to specify epitope sequence to search for!"
