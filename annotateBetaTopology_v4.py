import Tkinter
import tkFileDialog
import tkMessageBox
from pymol import cmd
import sys, urllib2, zlib
import os
import glob
import re
from collections import defaultdict
import dssp
 
def __init__(self):
	self.menuBar.addmenuitem('Plugin', 'command',
									 'Beta Topology Annotation v4',
									 label = 'Beta Topology Annotation v4',
									 command = lambda s=self : fetchAnnotationApp(s))

def infinite_defaultdict():
	return defaultdict(infinite_defaultdict)
###END infinite_deraultdict()###

def runDssp(pdbFile,dsspBinaryPath):
	dsspFilename = "temp.dssp"
	if os.path.isfile(dsspFilename):
		os.popen("rm dsspFilename")
	cmd = dsspBinaryPath + " " + pdbFile + " " + dsspFilename
	cmdResult = os.popen(cmd)
	print cmdResult.read()
	return(dsspFilename)

def getDsspSecStructHash(dsspFilename):

	mydssp = dssp.Dssp(dsspFilename)
	mydssp.setDsspData()
	dsspData = mydssp.dsspData
	dssp2pdbResnum = mydssp.dssp2pdbResnum

	return(dsspData,dssp2pdbResnum)
	

def loadMol(pdbFile):
	cmd.delete("all")
	try:
		print 'Loading ' + pdbFile
		cmd.load(pdbFile)
	except:
		print "Unexpected error:", sys.exc_info()[0]
		tkMessageBox.showerror('File I/O Error','Error while opening: ' + pdbFile)

def showBetaSheets():
	cmd.hide("everything")
	cmd.select("s","ss s")
	cmd.show("cartoon","s")
	cmd.color("grey","s")

def highlightBetaSheets():
	cmd.hide("everything")
	cmd.select("a","all")
	cmd.show("ribbon","a")
	cmd.color("forest","a")

	cmd.select("s","ss s")
	cmd.show("cartoon","s")
	cmd.color("grey","s")

def highlightSecStruct():
	cmd.hide("everything")
	cmd.select("a","all")
	cmd.show("ribbon","a")
	cmd.color("forest","a")

	cmd.select("strand","ss s")
	cmd.show("cartoon","strand")
	cmd.color("grey","strand")

	cmd.select("helix","ss h")
	#cmd.show("cartoon","helix")
	cmd.color("red","helix")


def highlightResidue(chain,resNum,representation,color,secStruct):
	cmd.select("res","resi " + resNum + " and chain " + chain)
	cmd.alter("res","ss=\'"+secStruct+"\'")
	cmd.show(representation,"res")
	cmd.color(color,"res")

def highlightDsspSecStruct(pdbFile,dsspBinaryPath):
	(dsspData,dssp2pdbResnum) = getDsspSecStructHash(runDssp(pdbFile,dsspBinaryPath))

	cmd.hide("everything")
	cmd.select("a","all")
	cmd.alter("a","ss=\'L\'")
	cmd.show("cartoon","a")
	cmd.color("deepblue","a")

	for chain in dsspData.keys():
		for pdbResNum in dsspData[chain].keys():
			print "Reading chain:" + chain + ", residue:" + pdbResNum + ", secstr:" + dsspData[chain][pdbResNum]["SEC_STRUCT"]
			if dsspData[chain][pdbResNum]["SEC_STRUCT"] == 'E':
				highlightResidue(chain,pdbResNum,"cartoon","grey",'S')
			elif dsspData[chain][pdbResNum]["SEC_STRUCT"] == 'B':
				highlightResidue(chain,pdbResNum,"cartoon","grey90",'S')
			elif dsspData[chain][pdbResNum]["SEC_STRUCT"] == 'H':
				highlightResidue(chain,pdbResNum,"cartoon","tv_red",'H')
			elif dsspData[chain][pdbResNum]["SEC_STRUCT"] == 'G':
				highlightResidue(chain,pdbResNum,"cartoon","magenta",'H')
			elif dsspData[chain][pdbResNum]["SEC_STRUCT"] == 'I':
				highlightResidue(chain,pdbResNum,"cartoon","pink",'H')
			elif dsspData[chain][pdbResNum]["SEC_STRUCT"] == 'T':
				highlightResidue(chain,pdbResNum,"cartoon","tv_green",'L')
			elif dsspData[chain][pdbResNum]["SEC_STRUCT"] == 'S':
				highlightResidue(chain,pdbResNum,"cartoon","limon",'L')
			else:
				highlightResidue(chain,pdbResNum,"cartoon","forest",'L')

	dsspData.clear()



def loadBetaSheets(pdbFile,button):
	loadMol(pdbFile)
	highlightBetaSheets()
	button.configure(bg = "red")

def loadSecStruct(pdbFile,dsspBinaryPath,button):
	loadMol(pdbFile)
	highlightDsspSecStruct(pdbFile,dsspBinaryPath)
	button.configure(bg = "red")


def getAnnotationHash(outputFilePath):

	outputFile = open(outputFilePath,'r')
	annotationHash = infinite_defaultdict()
	firstLine = ""
	for line in outputFile.readlines():
		line = line.replace("\n","")
		if line != "":
			if line.startswith("#"):
#				#PDB_ID	TRUE_BARREL	PSEUDO_BARREL	ONLY_BARREL	ALPHA_BETA_TOPOLOGY	NO_BETA	DOUBT	COMMENTS
				firstLine = line
			else:
				line = line.split("\t")
				annotationHash[line[0]] = {'TRUE_BARREL':line[1],'PSEUDO_BARREL':line[2],'ONLY_BARREL':line[3],'ALPHA_BETA_TOPOLOGY':line[4],'NO_BETA':line[5],'DOUBT':line[6],'COMMENTS':line[7]}

	outputFile.close()

	return(annotationHash,firstLine)


def updateAnnotationFile(buttons,vars_trueBarrel,vars_pseudoBarrel,vars_onlyBarrel,vars_alphaBetaTopology,vars_noBeta,vars_doubt,vars_comments,outputFilePath):
###READ THE EXISTING FILE. THIS WILL BE USED TO UPDATE THE EXISTING RECORDS.
	[annotationHash,firstLine] = getAnnotationHash(outputFilePath)

	for i in range(0,len(buttons)):
		if buttons[i].cget('bg') == "red":
			matchObj = re.match('\d+\.\s(.+)',buttons[i]['text'],re.S)
			pdbId = matchObj.group(1)
			print  "Updating annotationHash for " + pdbId
			annotationHash[pdbId] = {'TRUE_BARREL':str(vars_trueBarrel[i].get()),'PSEUDO_BARREL':str(vars_pseudoBarrel[i].get()),'ONLY_BARREL':str(vars_onlyBarrel[i].get()),'ALPHA_BETA_TOPOLOGY':vars_alphaBetaTopology[i].get(),'NO_BETA':str(vars_noBeta[i].get()),'DOUBT':str(vars_doubt[i].get()),'COMMENTS':vars_comments[i].get()}

	outputFile = open(outputFilePath,'w')
	outputFile.write(firstLine + "\n")
	for pdbId in sorted(annotationHash.keys()):
		outputFile.write(pdbId + "\t" + annotationHash[pdbId]['TRUE_BARREL'] + "\t" + annotationHash[pdbId]['PSEUDO_BARREL'] + "\t" + annotationHash[pdbId]['ONLY_BARREL'] + "\t" + annotationHash[pdbId]['ALPHA_BETA_TOPOLOGY'] + "\t" + annotationHash[pdbId]['NO_BETA'] + "\t" + annotationHash[pdbId]['DOUBT'] + "\t" + annotationHash[pdbId]['COMMENTS'] + "\n")

	outputFile.close()
	print "File  '" + outputFilePath + "'  updated successfully.\n"

def readConfig(configFilePath):
	configFile = open(configFilePath,'r')
	config = defaultdict()
	for line in configFile.readlines():
		line = line.replace("\n","")
		if line != "":
			matchObj = re.match('(.+)\t(.+)',line,re.S)
			config[matchObj.group(1)] = matchObj.group(2)

	return(config)

def fetchAnnotationApp(app):

###ASK INPUT DIRECTORY (OPTIONAL)
#	dirname = tkFileDialog.askdirectory(parent=app.root,initialdir="/home/harsh",title='Please select a directory')
#	if len(dirname ) > 0:
#		print "You chose %s" % dirname

###ASK FOR CONFIGURATION FILE
	configFilename = tkFileDialog.askopenfilename(title = "Select configuration file")
	print "Using configuration file: " + configFilename

###READ CONFIG FILE
	config = readConfig(configFilename)

	inputdir = config['inputdir']
	outputFilePath = config['outputFilePath']
	fileNumBegin = int(config['fileNumBegin'])
	if config['fileNumEnd'] == 'max':
		fileNumEnd = len(glob.glob( os.path.join(inputdir, '*.pdb')))
	else:
		fileNumEnd = int(config['fileNumEnd'])

	pdbFiles = sorted(glob.glob( os.path.join(inputdir, '*.pdb')))
	print "Displaying " + str(fileNumEnd-fileNumBegin+1) + " pdb files out of " + str(len(pdbFiles)) + " read."

	sizex = 680
	sizey = 494
	posx  = 700
	posy  = 100

	top = Tkinter.Tk()
	top.wm_geometry("%dx%d+%d+%d" % (sizex, sizey, posx, posy))
	top.wm_title("Beta Topology Annotation v3")

	headerFrame = Tkinter.Frame(top,width=660)
	headerFrame.place(x=10,y=10)

	Tkinter.Label(headerFrame, text = "     PDB ID    \n",bg="#AAAAAA",width=18).grid(row=2,column=1)
	Tkinter.Label(headerFrame, text = "  True  \n barrel ",bg="#CCCCCC").grid(row=2,column=2)
	Tkinter.Label(headerFrame, text = "Pseudo \n barrel ",bg="#AAAAAA").grid(row=2,column=3)
	Tkinter.Label(headerFrame, text = "Only \n barrel ",bg="#CCCCCC").grid(row=2,column=4)
	Tkinter.Label(headerFrame, text = "AB \n  topology  ",bg="#AAAAAA").grid(row=2,column=5)
	Tkinter.Label(headerFrame, text = "No Beta \n",bg="#CCCCCC").grid(row=2,column=6)
	Tkinter.Label(headerFrame, text = " Doubt \n",bg="#AAAAAA").grid(row=2,column=7)
	Tkinter.Label(headerFrame, text = "            Comments              \n",bg="#CCCCCC").grid(row=2,column=8)

	myframe = Tkinter.Frame(top,relief='groove',width=300,height=420,bd=1)
	myframe.place(x=10,y=45)

	canvas = Tkinter.Canvas(myframe)
	frame = Tkinter.Frame(canvas)
	myscrollbar=Tkinter.Scrollbar(myframe,orient="vertical",command=canvas.yview)
	canvas.configure(yscrollcommand=myscrollbar.set)

	myscrollbar.pack(side="right",fill="y")
	canvas.pack(side="left")
	canvas.create_window((0,0),window=frame,anchor='nw')
	frame.bind("<Configure>",lambda event: {canvas.configure(scrollregion=canvas.bbox("all"),width=640,height=400)})

###CREATE BUTTONS' ARRAY
	vars_trueBarrel = []
	vars_pseudoBarrel = []
	vars_onlyBarrel = []
	vars_alphaBetaTopology = []
	vars_noBeta = []
	vars_doubt = []
	vars_comments = []
	for i in range(0,fileNumEnd-fileNumBegin+1):
		vars_trueBarrel.append(Tkinter.IntVar(frame))
		vars_pseudoBarrel.append(Tkinter.IntVar(frame))
		vars_onlyBarrel.append(Tkinter.IntVar(frame))
		vars_alphaBetaTopology.append(Tkinter.StringVar(frame))
		vars_noBeta.append(Tkinter.IntVar(frame))
		vars_doubt.append(Tkinter.IntVar(frame))
		vars_comments.append(Tkinter.StringVar(frame))

	buttons = []
	checkButtons_trueBarrel = []
	checkButtons_pseudoBarrel = []
	checkButtons_onlyBarrel = []
	entries_alphaBetaTopology = []
	checkButtons_noBeta = []
	checkButtons_doubt = []
	entries_comments = []
	for i in range(0,fileNumEnd-fileNumBegin+1):
		buttons.append(Tkinter.Button(frame))
		checkButtons_trueBarrel.append(Tkinter.Checkbutton(frame, variable=vars_trueBarrel[i]))
		checkButtons_pseudoBarrel.append(Tkinter.Checkbutton(frame, variable=vars_pseudoBarrel[i]))
		checkButtons_onlyBarrel.append(Tkinter.Checkbutton(frame, variable=vars_onlyBarrel[i]))
		entries_alphaBetaTopology.append(Tkinter.Entry(frame, textvariable=vars_alphaBetaTopology[i], width=8))
		checkButtons_noBeta.append(Tkinter.Checkbutton(frame, variable=vars_noBeta[i]))
		checkButtons_doubt.append(Tkinter.Checkbutton(frame, variable=vars_doubt[i]))
		entries_comments.append(Tkinter.Entry(frame, textvariable=vars_comments[i], width=20))

	goButton = Tkinter.Button(top, text = "Go", command = lambda:updateAnnotationFile(buttons,vars_trueBarrel,vars_pseudoBarrel,vars_onlyBarrel,vars_alphaBetaTopology,vars_noBeta,vars_doubt,vars_comments,config['outputFilePath']))

###SET BUTTON PARAMETERS
	[annotationHash,firstLine] = getAnnotationHash(outputFilePath)

	cnt = 0
	for i in range(fileNumBegin-1,fileNumEnd):
		matchObj = re.match('.*/(.+)\.pdb$',pdbFiles[i],re.S)
		pdbId = matchObj.group(1)

		buttons[cnt].configure(text = str(i+1) + ". " + pdbId, width=15)
		buttons[cnt].configure(command = lambda j=i,cnt=cnt: loadSecStruct(pdbFiles[j],config['dsspBinaryPath'],buttons[cnt]))

		if pdbId in annotationHash.keys():
			buttons[cnt].configure(bg='green')
			if annotationHash[pdbId]['TRUE_BARREL'] == '1':
				checkButtons_trueBarrel[cnt].select()
			if annotationHash[pdbId]['PSEUDO_BARREL'] == '1':
				checkButtons_pseudoBarrel[cnt].select()
			if annotationHash[pdbId]['ONLY_BARREL'] == '1':
				checkButtons_onlyBarrel[cnt].select()
			vars_alphaBetaTopology[cnt].set(annotationHash[pdbId]['ALPHA_BETA_TOPOLOGY'])
			if annotationHash[pdbId]['NO_BETA'] == '1':
				checkButtons_noBeta[cnt].select()
			if annotationHash[pdbId]['DOUBT'] == '1':
				checkButtons_doubt[cnt].select()
			vars_comments[cnt].set(annotationHash[pdbId]['COMMENTS'])

		cnt += 1

###RENDER BUTTONS
	Tkinter.Label(frame, text = "      ",width=14).grid(row=0,column=1)
	Tkinter.Label(frame, text = "    ",width=5).grid(row=0,column=2)
	Tkinter.Label(frame, text = "      ",width=6).grid(row=0,column=3)
	Tkinter.Label(frame, text = "      ",width=6).grid(row=0,column=4)
	Tkinter.Label(frame, text = "    ",width=5).grid(row=0,column=5)
	Tkinter.Label(frame, text = "       ",width=6).grid(row=0,column=6)
	Tkinter.Label(frame, text = "     ",width=6).grid(row=0,column=7)

	for i in range(1,len(buttons)+1):
		buttons[i-1].grid(row=i, column=1)
		checkButtons_trueBarrel[i-1].grid(row=i,column=2)
		checkButtons_pseudoBarrel[i-1].grid(row=i,column=3)
		checkButtons_onlyBarrel[i-1].grid(row=i,column=4)
		entries_alphaBetaTopology[i-1].grid(row=i,column=5)
		checkButtons_noBeta[i-1].grid(row=i,column=6)
		checkButtons_doubt[i-1].grid(row=i,column=7)
		entries_comments[i-1].grid(row=i,column=8)


	goButton.place(x=628,y=458)

	top.mainloop()

cmd.extend('loadSecStruct', loadSecStruct)

