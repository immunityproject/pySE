import sys
import uuid
import os


DEBUG = False


class amino:
	location = 0
	code = "X"
	chain = "X"
	c_loc = 0


	def __init__(self, location, code, chain, c_loc):
		self.location = location
		self.code = code
		self.chain = chain
		self.c_loc = c_loc


codes = {"ALA": "A",
			"ARG": "R",
			"ASN": "N",
			"ASP": "D",
			"ASX": "B",
			"CYS": "C",
			"GLU": "E",
			"GLN": "Q",
			"GLX": "Z",
			"GLY": "G",
			"HIS": "H",
			"ILE": "I",
			"LEU": "L",
			"LYS": "K",
			"MET": "M",
			"PHE": "F",
			"PRO": "P",
			"SER": "S",
			"THR": "T",
			"TRP": "W",
			"TYR": "Y",
			"VAL": "V",
			"TPO": "X"}


def readPDB_seqres(filename):
	amino_dic = {}
	f = open(filename, 'r')

	location = 0
	code = ""
	chain = ""
	c_loc = 0

	for line in f:
		line_array = line.split()
		if len(line_array) > 0:
			if line_array[0] == "SEQRES":
				if chain != line_array[2]:
					c_loc = 0
				chain = line_array[2]
				for i in (range(4, len(line_array))):
					location = location + 1
					code = codes.get(line_array[i])
					c_loc = c_loc + 1
					amino_dic[location] = amino(location, code, chain, c_loc)

	return amino_dic


def readPDB_atom(filename):
	amino_dic = {}
	f = open(filename, 'r')

	location = -1
	code = ""
	chain = ""
	c_loc = 0

	first_loc = -1
	offset = 0

	last_c_loc = 0

	for line in f:
		line_array = line.split()
		if len(line_array) > 0:
			if line_array[0] == "ATOM":

				# skip lines for atom repetition.  We don't care about their other details.
				# This is degenerate in 1 atom length chains!
				try:
					if c_loc != int(line_array[5]):
						if chain != line_array[4]:
							location = location + int(line_array[5]) - 1
							last_c_loc = location - 1

							# is this valid? Chains start at non-1
							if (first_loc == -1):
								first_loc = int(line_array[5])
								# undo location and last_c_loc for first amino.
								location = first_loc - 1
								last_c_loc = location

						location = location + 1
						code = codes.get(line_array[3])
						chain = line_array[4]
						c_loc = int(line_array[5])
						diff = c_loc - last_c_loc

						if (DEBUG):
							print "DIFF IS:" + str(diff)

						# fill in missing amino data with "X".
						if (diff > 1):
							for y in range(last_c_loc + 1, c_loc):
								amino_dic[location - 1] = amino(location, "X", chain, y)
								location = location + 1
						last_c_loc = c_loc
						try:
							if (DEBUG):
								print "location:" + str(location),
								print ",code:" + code,
								print ",chain:" + chain,
								print ",c_loc:" + str(c_loc),
								print",first_loc:" + str(first_loc),

							amino_dic[location - 1] = amino(location, code, chain, c_loc)

							if (DEBUG):
								print " | rechecking:",
								print amino_dic[location - 1].code
						except Exception, e:
							print "Found an error on the following line:"
							print line

						# else:
						# print "I skipped a line because it didn't start with ATOM"
				except Exception:
					print "The following line was malformed:"
					print line
	return first_loc, amino_dic


def generate_foldx_posscan(amino_dic, start, stop):
	print "GENERATING POSITIONSCAN CONFIGURATION STRING"
	return_string = ""
	amino_string = ""

	for i in range(start - 1, stop):
		amino_string = amino_string + amino_dic[i].code
		return_string = return_string + amino_dic[i].code + amino_dic[i].chain + str(amino_dic[i].c_loc) + "a,"

	print amino_string
	print return_string[0:len(return_string) - 1]


def generate_foldx_buildmodel_single(amino_dic, loc, chains, code):
	line = ""
	for k in range(0, chains):
		if (len(line) > 0):
			line = line + ","
		line = line + amino_dic[loc].code + chr(65 + k) + str(loc + 1) + code
	line = line + ";\n"
	return line


def generate_foldx_buildmodel(amino_dic, start, stop, chains):
	count = 0

	return_line = ""

	for i in range(start - 1, stop):
		for code in codes.values():
			line = generate_foldx_buildmodel_single(amino_dic, i, chains, code)
			return_line = return_line + line
			count = count + 1

	return count, return_line


def generate_foldx_buildmodel_file(amino_dic, start, stop, chains):
	print "GENERATING POSITIONSCAN CONFIGURATION STRING"
	f = open('./individual_list.txt', 'w+')

	count, return_line = generate_foldx_buildmodel(amino_dic, start, stop, chains)

	f.write(return_line)
	print "Wrote " + str(count) + " lines to individual_list.txt"
	f.close()


def generate_foldx_buildmodel_dirs(amino_dic, start, stop, chains, pdb, root):
	for i in range(start - 1, stop):
		for code in codes.values():
			line = generate_foldx_buildmodel_single(amino_dic, i, chains, code)

			dir_name = "foldxbm-" + str(uuid.uuid4())
			full_loc = root + "/" + dir_name
			os.makedirs(full_loc)
			f = open(full_loc + '/individual_list.txt', 'w+')
			f2 = open(full_loc + '/list.txt', 'w+')
			f2.write(os.path.basename(pdb) + "\n")
			f2.close()
			f.write(line)
			f.close()


# This is broken needs to increment over dict values like below TODO
def dump_chain(first_loc, amino_dic):
	full_chain = ""
	print "LENGTH:" + str(len(amino_dic))
	for i in range(len(amino_dic)):
		try:
			if amino_dic[first_loc + i - 1].chain == "A":
				full_chain = full_chain + amino_dic[first_loc + i - 1].code
		except KeyError:
			pass

	print full_chain


def dump_all_chain(first_loc, amino_dic):
	full_chain = ""
	hits = 0
	current_try = 0

	print "LENGTH:" + str(len(amino_dic))
	while hits < len(amino_dic):
		try:
			if (DEBUG):
				print current_try,
				print ":",
				print amino_dic[first_loc + current_try - 1].code
			# print "Trying position:" + str(first_loc+current_try-1)
			if (amino_dic[first_loc + current_try - 1].code):
				hits = hits + 1
				full_chain = full_chain + amino_dic[first_loc + current_try - 1].code

			current_try = current_try + 1

		except KeyError:
			print "keyerror"
			current_try = current_try + 1
			pass

	print full_chain


def main():
	#	for i in xrange(15):
	#		print uuid.uuid4()
	#
	root = "./"
	generateDistributableDirs = 0

	if len(sys.argv) < 4:
		print "Improper arguments.  Usage: " + sys.argv[
			0] + " <PDB> <start position> <stop position> <isBuildModel? 0 or 1> [NumberChains (0 or 1)] [generate distributable dirs (0 or 1)] [dir root]"
		exit();
	if int(sys.argv[2]) <= 0:
		print "Start must be greater than zero.  Indices start at 1, not zero."
		exit();

	if int(sys.argv[4]):
		if int(sys.argv[5]) <= 0:
			print "Chains is always greater than or equal to 1."
			exit();
		else:
			chains = int(sys.argv[5])

	if len(sys.argv) > 6:
		generateDistributableDirs = int(sys.argv[6])

	if len(sys.argv) > 7:
		root = sys.argv[7]
		if not os.path.exists(root):
			os.makedirs(root)

	pdb_name = sys.argv[1]
	start = int(sys.argv[2])
	stop = int(sys.argv[3])
	isBuildModel = int(sys.argv[4])

	# amino_dic = readPDB_seqres(pdb_name)
	first_loc, amino_dic = readPDB_atom(pdb_name)

	print "dumping chain..."
	print "\n"
	print dump_all_chain(first_loc, amino_dic)
	print "\n"

	if (not isBuildModel):
		generate_foldx_posscan(amino_dic, start, stop)
	else:
		if (not generateDistributableDirs):
			generate_foldx_buildmodel_file(amino_dic, start, stop, chains)
		else:
			generate_foldx_buildmodel_dirs(amino_dic, start, stop, chains, pdb_name, root)


if __name__ == '__main__':
	main();
