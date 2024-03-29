import re
from collections import defaultdict

def infinite_defaultdict():
	return defaultdict(infinite_defaultdict)
###END infinite_deraultdict()###


class Dssp :

	def __init__(self,filename):
		self.filename = filename

###DEFINE GLOBAL VARIABLES FOR CLASS Dssp
		self.dsspData = infinite_defaultdict()
		self.dssp2pdbResnum = infinite_defaultdict()
		self.aaCode1to3 = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}

###END __init__()###

	def setDsspData(self):

		print "Class Dssp :: Opening file",self.filename
		dsspfile = open(self.filename,'r')
		lines =  dsspfile.readlines()
		
		beginData = 0
		for line in lines:
			#print line

			if beginData == 1:
				"""
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA 
   12   12 A H  H  < S+     0   0   12     -4,-1.9   160,-3.6   159,-0.1     2,-0.5   0.501  96.4  91.0-140.8  -8.6    7.2   -6.4   15.0
   13   13 A M  B  <  +a  172   0A  13     -4,-1.8   160,-0.2    -5,-0.4    -1,-0.1  -0.858  24.4 147.9-104.5 129.0    7.9   -3.2   16.7
  210   85 B R  E     -BC 168 223A 120     13,-2.1    13,-3.2    -2,-0.6   -42,-0.2  -0.964  33.3-111.8-120.3 132.7   -3.0  -15.2   15.6
  211   86 B E  E     - C   0 222A  65    -44,-2.4    11,-0.3    -2,-0.5     2,-0.1  -0.293  33.7-131.0 -55.2 139.7   -5.5  -15.7   12.8
   14   14 A D        -     0   0   33    158,-2.3     2,-0.1    -2,-0.5   161,-0.1  -0.542  29.5-166.9-148.8  77.6    5.2   -1.0   18.3
   24   32BA S              0   0  154      1,-0.2     2,-0.3     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0  22.6    8.4   16.3   32.3
				"""
#==============================   13==   13 ==A= =M=  =B==  <  +a == 172==   0==A==  13= =    -4,-1.8==   160,-0.2==    -5,-0.4==    -1,-0.1=  =-0.858==  24.4== 147.9==-104.5== 129.0=    7.9   -3.2   16.7
#==============================   14==   14A==A= =D=  = ==     -  ==   0==   0== ==  33= =   158,-2.3==     2,-0.1==    -2,-0.5==   161,-0.1=  =-0.542==  29.5==-166.9==-148.8==  77.6=    5.2   -1.0   18.3
#=============================(.....)(......)(.).(.)..(.)(........)(....)(....)(.)(....).(...........)(...........)(...........)(...........)..(......)(......)(......)(......)(......).*
				matchObj = re.match(r'(.....)(......)(.).(.)..(.)(........)(....)(....)(.)(....).(...........)(...........)(...........)(...........)..(......)(......)(......)(......)(......).*\n',line,re.S)
				if matchObj.group(4) != "!":
					DSSP_RES_NUM = int(matchObj.group(1).strip())
					PDB_RES_NUM = matchObj.group(2).strip()
					CHAIN = matchObj.group(3).strip()
					if matchObj.group(4).strip() in self.aaCode1to3:
						RES_NAME = self.aaCode1to3[matchObj.group(4).strip()]
					else:
						RES_NAME = 'UNK'
					SEC_STRUCT = matchObj.group(5)
					SEC_STRUCT_INFO = matchObj.group(6)
					BP1 = int(matchObj.group(7).strip())
					BP2 = int(matchObj.group(8).strip())
					SHEET_ID = matchObj.group(9).strip()
					ACC = float(matchObj.group(10).strip())
					NH_O_1 = matchObj.group(11).strip()
					O_HN_1 = matchObj.group(12).strip()
					NH_O_2 = matchObj.group(13).strip()
					O_HN_2 = matchObj.group(14).strip()
					TCO = float(matchObj.group(15).strip())
					KAPPA = float(matchObj.group(16).strip())
					ALPHA = float(matchObj.group(17).strip())
					PHI = float(matchObj.group(18).strip())
					PSI = float(matchObj.group(19).strip())


					self.dsspData[CHAIN][PDB_RES_NUM] =\
					{
						'DSSP_RES_NUM':DSSP_RES_NUM,
						'RES_NAME':RES_NAME,
						'SEC_STRUCT':SEC_STRUCT,
						'SEC_STRUCT_INFO':SEC_STRUCT_INFO,
						'BP1':BP1,
						'BP2':BP2,
						'SHEET_ID':SHEET_ID,
						'ACC':ACC,
						'NH_O_1':NH_O_1,
						'O_HN_1':O_HN_1,
						'NH_O_2':NH_O_2,
						'O_HN_2':O_HN_2,
						'TCO':TCO,
						'KAPPA':KAPPA,
						'ALPHA':ALPHA,
						'PHI':PHI,
						'PSI':PSI
					}

					self.dssp2pdbResnum[DSSP_RES_NUM] = PDB_RES_NUM

			if line.startswith("  #"):
				beginData = 1

		dsspfile.close()





