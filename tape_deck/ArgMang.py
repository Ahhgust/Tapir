import io
import os
import sys
import struct
import subprocess

try:
	import tkinter
	import tkinter.filedialog
	import tkinter.messagebox
	import tkinter.scrolledtext
except ImportError as e:
	sys.stderr.write("Need a version of python with TKinter.")
	raise

import whodunargs


class ProgSetInfo:
	'''Information on a program set.'''
	def __init__(self,fromStr):
		textL = struct.unpack(">Q", fromStr.read(8))[0]
		self.name = str(fromStr.read(textL), "utf-8")
		'''The name of the set.'''
		textL = struct.unpack(">Q", fromStr.read(8))[0]
		self.summary = str(fromStr.read(textL), "utf-8")
		'''A summary of the set.'''
		numProg = struct.unpack(">Q", fromStr.read(8))[0]
		self.progNames = []
		'''The names of the programs'''
		self.progSummaries = []
		'''Quick summaries of the programs'''
		for i in range(numProg):
			textL = struct.unpack(">Q", fromStr.read(8))[0]
			tmpName = str(fromStr.read(textL), "utf-8")
			textL = struct.unpack(">Q", fromStr.read(8))[0]
			tmpSummary = str(fromStr.read(textL), "utf-8")
			self.progNames.append(tmpName)
			self.progSummaries.append(tmpSummary)


class ArgInfo:
	'''Information on a command line argument.'''
	def __init__(self,fromStr):
		textL = struct.unpack(">Q", fromStr.read(8))[0]
		self.name = str(fromStr.read(textL), "utf-8")
		'''The name of the argument.'''
		sigils = []
		numSig = struct.unpack(">Q", fromStr.read(8))[0]
		for i in range(numSig):
			textL = struct.unpack(">Q", fromStr.read(8))[0]
			sigils.append(str(fromStr.read(textL), "utf-8"))
		self.sigils = sigils
		'''The markers for this argument.'''
		textL = struct.unpack(">Q", fromStr.read(8))[0]
		self.summary = str(fromStr.read(textL), "utf-8")
		'''A summary of the argument.'''
		textL = struct.unpack(">Q", fromStr.read(8))[0]
		self.description = str(fromStr.read(textL), "utf-8")
		'''A longform description of the argument.'''
		textL = struct.unpack(">Q", fromStr.read(8))[0]
		self.usage = str(fromStr.read(textL), "utf-8")
		'''An example usage of the argument.'''
		self.isPublic = fromStr.read(1)[0] != 0
		'''Whether this argument should be visible.'''
		textL = struct.unpack(">Q", fromStr.read(8))[0]
		self.mainFlavor = str(fromStr.read(textL), "utf-8")
		'''The type of argument this is.'''
		textL = struct.unpack(">Q", fromStr.read(8))[0]
		self.subFlavor = str(fromStr.read(textL), "utf-8")
		'''The sub-type.'''
		extL = struct.unpack(">Q", fromStr.read(8))[0]
		self.extras = fromStr.read(extL)
		'''The extra data.'''


class ProgramInfo:
	'''Information on a program.'''
	def __init__(self,fromStr):
		textL = struct.unpack(">Q", fromStr.read(8))[0]
		self.name = str(fromStr.read(textL), "utf-8")
		'''The name of the set.'''
		textL = struct.unpack(">Q", fromStr.read(8))[0]
		self.summary = str(fromStr.read(textL), "utf-8")
		'''A summary of the set.'''
		textL = struct.unpack(">Q", fromStr.read(8))[0]
		self.description = str(fromStr.read(textL), "utf-8")
		'''A longform description of the program.'''
		textL = struct.unpack(">Q", fromStr.read(8))[0]
		self.usage = str(fromStr.read(textL), "utf-8")
		'''An example usage of the program.'''
		numArg = struct.unpack(">Q", fromStr.read(8))[0]
		theArgs = []
		for i in range(numArg):
			theArgs.append(ArgInfo(fromStr))
		self.arguments = theArgs
		'''The arguments this takes.'''


def getProgramInfos(progPath,subPName = None):
	'''
	Get info for a program.
	@param progPath: The path to the program, and any prefix arguments.
	@param subPName: The name of the subprogram of interest, or None if not a set.
	@return: The info.
	'''
	runArgs = progPath[:]
	if not (subPName is None):
		runArgs.append(subPName)
	runArgs.append("--help_argdump")
	curRunP = subprocess.Popen(runArgs,stdout=subprocess.PIPE,stdin=subprocess.DEVNULL)
	allArgs = ProgramInfo(curRunP.stdout)
	curRunP.stdout.close()
	if curRunP.wait() != 0:
		raise IOError("Problem getting argument info.")
	return allArgs


def getProgramSetInfos(progPath):
	'''
	Get info for a program set.
	@param progPath: The path to the program.
	@param subPName: The name of the subprogram of interest, or None if not a set.
	@return: The info.
	'''
	runArgs = progPath[:]
	runArgs.append("--help_argdump")
	curRunP = subprocess.Popen(runArgs,stdout=subprocess.PIPE,stdin=subprocess.DEVNULL)
	allArgs = ProgSetInfo(curRunP.stdout)
	curRunP.stdout.close()
	if curRunP.wait() != 0:
		raise IOError("Problem getting set info.")
	return allArgs


class ArgumentFlavor:
	def manpageSynopMang(self,argD):
		'''
		Make some text for an argument synopsis in a manpage.
		@param argD: The argument data.
		'''
		raise ValueError("Not implemented for argument type")
	def makeArgGUI(self,forGui,argD):
		'''
		Make a gui for an argument.
		@param forGui: The gui to add to.
		@param argD: The argument data.
		'''
		raise ValueError("Not implemented for argument type")
	def getGUIArgTexts(self,forGui,forData):
		'''
		Get the argument texts for a given set of arguments.
		@param forGui: The gui this is for.
		@param forData: The data this made in makeArgGUI.
		@return: A list of argument
		'''
		raise ValueError("Not implemented for argument type")


class ArgumentFlavorFlag(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "\\fB" + manSan(argD.sigils[0]) + "\\fR"
	def makeArgGUI(self,forGui,argD):
		theV = tkinter.IntVar()
		theCB = tkinter.Checkbutton(forGui.myCanvas, text=argD.name, variable=theV)
		theCB.grid(column = 0, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theV, theCB, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		if forData[1].get() != 0:
			toRet.append(forData[0].sigils[0])
		return toRet


class ArgumentFlavorEnum(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "\\fB" + manSan(argD.sigils[0]) + "\\fR"
	def makeArgGUI(self,forGui,argD):
		# see if the class is in the handled classes
		if not ("enum_class" in forGui.valueCache):
			forGui.valueCache["enum_class"] = {}
		classNameL = struct.unpack(">Q", argD.extras[0:8])[0]
		className = str(argD.extras[8:(classNameL+8)], "utf-8")
		if className in forGui.valueCache["enum_class"]:
			tgtStrVar = forGui.valueCache["enum_class"][className]
		else:
			tgtStrVar = tkinter.StringVar(forGui.master, argD.name)
			forGui.valueCache["enum_class"][className] = tgtStrVar
		# make the radio button
		theCB = tkinter.Radiobutton(forGui.myCanvas, text=argD.name, variable=tgtStrVar, value=argD.name)
		theCB.grid(column = 0, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, tgtStrVar, theCB, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		if forData[1].get() == forData[0].name:
			toRet.append(forData[0].sigils[0])
		return toRet


class ArgumentFlavorInt(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "\\fB" + manSan(argD.sigils[0]) + "\\fR \\fI###\\fR"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theBox = tkinter.Entry(forGui.myCanvas)
		theBox.grid(column = 1, row = forGui.gridR)
		theBox.delete(0, tkinter.END)
		initValue = struct.unpack(">q", argD.extras[0:8])[0]
		theBox.insert(0, repr(initValue))
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curVal = forData[2].get().strip()
		if len(curVal) > 0:
			toRet.append(forData[0].sigils[0])
			toRet.append(forData[2].get().strip())
		return toRet


class ArgumentFlavorIntVec(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "[\\fB" + manSan(argD.sigils[0]) + "\\fR \\fI###\\fR]*"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		theBox = tkinter.scrolledtext.ScrolledText(forGui.myCanvas, height = 4, width = 40)
		theBox.grid(column = 0, row = forGui.gridR, columnspan = 2, rowspan = 4)
		forGui.gridR = forGui.gridR + 4
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curRetOpts = forData[2].get('1.0', tkinter.END).split("\n")
		for cv in curRetOpts:
			if len(cv.strip()) == 0:
				continue
			toRet.append(forData[0].sigils[0])
			toRet.append(cv.strip())
		return toRet


class ArgumentFlavorIntGreedVec(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "[\\fB" + manSan(argD.sigils[0]) + "\\fR \\fI###*\\fR]"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		theBox = tkinter.scrolledtext.ScrolledText(forGui.myCanvas, height = 4, width = 40)
		theBox.grid(column = 0, row = forGui.gridR, columnspan = 2, rowspan = 4)
		forGui.gridR = forGui.gridR + 4
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curRetOpts = forData[2].get('1.0', tkinter.END).split("\n")
		toRet.append(forData[0].sigils[0])
		for cv in curRetOpts:
			if len(cv.strip()) == 0:
				continue
			toRet.append(cv.strip())
		return toRet


class ArgumentFlavorFloat(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "\\fB" + manSan(argD.sigils[0]) + "\\fR \\fI###.###\\fR"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theBox = tkinter.Entry(forGui.myCanvas)
		theBox.grid(column = 1, row = forGui.gridR)
		theBox.delete(0, tkinter.END)
		initValue = struct.unpack(">d", argD.extras[0:8])[0]
		theBox.insert(0, repr(initValue))
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curVal = forData[2].get().strip()
		if len(curVal) > 0:
			toRet.append(forData[0].sigils[0])
			toRet.append(forData[2].get().strip())
		return toRet


class ArgumentFlavorFloatVec(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "[\\fB" + manSan(argD.sigils[0]) + "\\fR \\fI###.###\\fR]*"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		theBox = tkinter.scrolledtext.ScrolledText(forGui.myCanvas, height = 4, width = 40)
		theBox.grid(column = 0, row = forGui.gridR, columnspan = 2, rowspan = 4)
		forGui.gridR = forGui.gridR + 4
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curRetOpts = forData[2].get('1.0', tkinter.END).split("\n")
		for cv in curRetOpts:
			if len(cv.strip()) == 0:
				continue
			toRet.append(forData[0].sigils[0])
			toRet.append(cv.strip())
		return toRet


class ArgumentFlavorFloatGreedVec(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "[\\fB" + manSan(argD.sigils[0]) + "\\fR \\fI###.###*\\fR]"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		theBox = tkinter.scrolledtext.ScrolledText(forGui.myCanvas, height = 4, width = 40)
		theBox.grid(column = 0, row = forGui.gridR, columnspan = 2, rowspan = 4)
		forGui.gridR = forGui.gridR + 4
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curRetOpts = forData[2].get('1.0', tkinter.END).split("\n")
		toRet.append(forData[0].sigils[0])
		for cv in curRetOpts:
			if len(cv.strip()) == 0:
				continue
			toRet.append(cv.strip())
		return toRet


class ArgumentFlavorString(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "\\fB" + manSan(argD.sigils[0]) + "\\fR \\fITEXT\\fR"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theBox = tkinter.Entry(forGui.myCanvas)
		theBox.grid(column = 1, row = forGui.gridR)
		theBox.delete(0, tkinter.END)
		initLen = struct.unpack(">Q", argD.extras[0:8])[0]
		initValue = str(argD.extras[8:(initLen+8)], "utf-8")
		theBox.insert(0, initValue)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curVal = forData[2].get().strip()
		if len(curVal) > 0:
			toRet.append(forData[0].sigils[0])
			toRet.append(curVal)
		return toRet


class ArgumentFlavorStringVec(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "[\\fB" + manSan(argD.sigils[0]) + "\\fR \\fITEXT\\fR]*"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		theBox = tkinter.scrolledtext.ScrolledText(forGui.myCanvas, height = 4, width = 40)
		theBox.grid(column = 0, row = forGui.gridR, columnspan = 2, rowspan = 4)
		forGui.gridR = forGui.gridR + 4
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curRetOpts = forData[2].get('1.0', tkinter.END).split("\n")
		for cv in curRetOpts:
			if len(cv.strip()) == 0:
				continue
			toRet.append(forData[0].sigils[0])
			toRet.append(cv.strip())
		return toRet


class ArgumentFlavorStringGreedVec(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "[\\fB" + manSan(argD.sigils[0]) + "\\fR \\fITEXT*\\fR]"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		theBox = tkinter.scrolledtext.ScrolledText(forGui.myCanvas, height = 4, width = 40)
		theBox.grid(column = 0, row = forGui.gridR, columnspan = 2, rowspan = 4)
		forGui.gridR = forGui.gridR + 4
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curRetOpts = forData[2].get('1.0', tkinter.END).split("\n")
		curSigOn = 0
		for cv in curRetOpts:
			cvs = cv.strip()
			if len(cvs) == 0:
				continue
			if curSigOn == 0:
				toRet.append(forData[0].sigils[0])
				curSigOn = 1
			if (curSigOn > 1) and (len(cvs) >= 2) and (cvs[0] == '-') and (cvs[1] == '-'):
				toRet.append("--.-")
				toRet.append(forData[0].sigils[0])
				toRet.append(cvs)
				curSigOn = 0
			else:
				toRet.append(cvs)
				curSigOn = 2
		return toRet


def parseFileExtensions(forExtra):
	numExts = struct.unpack(">Q",forExtra[0:8])[0]
	extBytes = forExtra[8:]
	allExts = []
	for i in range(numExts):
		curELen = struct.unpack(">Q", extBytes[0:8])[0]
		curETxt = str(extBytes[8:(8+curELen)],"utf-8")
		extBytes = extBytes[8+curELen:]
		allExts.append(curETxt)
	return (allExts, extBytes)


class ArgumentFlavorStringFileRead(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "\\fB" + manSan(argD.sigils[0]) + "\\fR \\fIFILE\\fR"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theBox = tkinter.Entry(forGui.myCanvas)
		theBox.grid(column = 1, row = forGui.gridR)
		theBox.delete(0, tkinter.END)
		initLen = struct.unpack(">Q", argD.extras[0:8])[0]
		initValue = str(argD.extras[8:(initLen+8)], "utf-8")
		theBox.insert(0, initValue)
		allExts = parseFileExtensions(argD.extras[initLen+8:])[0]
		theBtn = None
		if len(allExts) > 0:
			theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateSingleFileName(theBox, tkinter.filedialog.askopenfilename(title = argD.name, filetypes = tuple([(cv, "*" + cv) for cv in allExts]) )))
		else:
			theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateSingleFileName(theBox, tkinter.filedialog.askopenfilename(title = argD.name )))
		theBtn.grid(column = 2, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theBtn, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curVal = forData[2].get().strip()
		if len(curVal) > 0:
			toRet.append(forData[0].sigils[0])
			toRet.append(curVal)
		return toRet


class ArgumentFlavorStringVecFileRead(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "[\\fB" + manSan(argD.sigils[0]) + "\\fR \\fIFILE\\fR]*"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		theBox = tkinter.scrolledtext.ScrolledText(forGui.myCanvas, height = 4, width = 40)
		theBox.grid(column = 1, row = forGui.gridR, columnspan = 2, rowspan = 4)
		forGui.gridR = forGui.gridR + 4
		allExts = parseFileExtensions(argD.extras)[0]
		theBtn = None
		if len(allExts) > 0:
			theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateMultiFileName(theBox, tkinter.filedialog.askopenfilenames(title = argD.name, filetypes = tuple([(cv, "*" + cv) for cv in allExts]) )))
		else:
			theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateMultiFileName(theBox, tkinter.filedialog.askopenfilenames(title = argD.name )))
		theBtn.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theBtn, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curRetOpts = forData[2].get('1.0', tkinter.END).split("\n")
		for cv in curRetOpts:
			if len(cv.strip()) == 0:
				continue
			toRet.append(forData[0].sigils[0])
			toRet.append(cv.strip())
		return toRet


class ArgumentFlavorStringGreedVecFileRead(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "[\\fB" + manSan(argD.sigils[0]) + "\\fR \\fIFILE*\\fR]"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		theBox = tkinter.scrolledtext.ScrolledText(forGui.myCanvas, height = 4, width = 40)
		theBox.grid(column = 1, row = forGui.gridR, columnspan = 2, rowspan = 4)
		forGui.gridR = forGui.gridR + 4
		allExts = parseFileExtensions(argD.extras)[0]
		theBtn = None
		if len(allExts) > 0:
			theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateMultiFileName(theBox, tkinter.filedialog.askopenfilenames(title = argD.name, filetypes = tuple([(cv, "*" + cv) for cv in allExts]) )))
		else:
			theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateMultiFileName(theBox, tkinter.filedialog.askopenfilenames(title = argD.name )))
		theBtn.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theBtn, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curRetOpts = forData[2].get('1.0', tkinter.END).split("\n")
		curSigOn = 0
		for cv in curRetOpts:
			cvs = cv.strip()
			if len(cvs) == 0:
				continue
			if curSigOn == 0:
				toRet.append(forData[0].sigils[0])
				curSigOn = 1
			if (curSigOn > 1) and (len(cvs) >= 2) and (cvs[0] == '-') and (cvs[1] == '-'):
				toRet.append("--.-")
				toRet.append(forData[0].sigils[0])
				toRet.append(cvs)
				curSigOn = 0
			else:
				toRet.append(cvs)
				curSigOn = 2
		return toRet


class ArgumentFlavorStringFileWrite(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "\\fB" + manSan(argD.sigils[0]) + "\\fR \\fIFILE\\fR"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theBox = tkinter.Entry(forGui.myCanvas)
		theBox.grid(column = 1, row = forGui.gridR)
		theBox.delete(0, tkinter.END)
		initLen = struct.unpack(">Q", argD.extras[0:8])[0]
		initValue = str(argD.extras[8:(initLen+8)], "utf-8")
		theBox.insert(0, initValue)
		allExts = parseFileExtensions(argD.extras[initLen+8:])[0]
		theBtn = None
		if len(allExts) > 0:
			theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateSingleFileName(theBox, tkinter.filedialog.asksaveasfilename(title = argD.name, filetypes = tuple([(cv, "*" + cv) for cv in allExts]), defaultextension="*.*" )))
		else:
			theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateSingleFileName(theBox, tkinter.filedialog.asksaveasfilename(title = argD.name )))
		theBtn.grid(column = 2, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theBtn, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curVal = forData[2].get().strip()
		if len(curVal) > 0:
			toRet.append(forData[0].sigils[0])
			toRet.append(curVal)
		return toRet


class ArgumentFlavorStringVecFileWrite(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "[\\fB" + manSan(argD.sigils[0]) + "\\fR \\fIFILE\\fR]*"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		theBox = tkinter.scrolledtext.ScrolledText(forGui.myCanvas, height = 4, width = 40)
		theBox.grid(column = 1, row = forGui.gridR, columnspan = 2, rowspan = 4)
		forGui.gridR = forGui.gridR + 4
		allExts = parseFileExtensions(argD.extras)[0]
		theBtn = None
		if len(allExts) > 0:
			theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateMultiFileName(theBox, tkinter.filedialog.asksaveasfilename(title = argD.name, filetypes = tuple([(cv, "*" + cv) for cv in allExts]), defaultextension="*.*" )))
		else:
			theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateMultiFileName(theBox, tkinter.filedialog.asksaveasfilename(title = argD.name )))
		theBtn.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theBtn, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curRetOpts = forData[2].get('1.0', tkinter.END).split("\n")
		for cv in curRetOpts:
			if len(cv.strip()) == 0:
				continue
			toRet.append(forData[0].sigils[0])
			toRet.append(cv.strip())
		return toRet


class ArgumentFlavorStringGreedVecFileWrite(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "[\\fB" + manSan(argD.sigils[0]) + "\\fR \\fIFILE*\\fR]"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		theBox = tkinter.scrolledtext.ScrolledText(forGui.myCanvas, height = 4, width = 40)
		theBox.grid(column = 1, row = forGui.gridR, columnspan = 2, rowspan = 4)
		forGui.gridR = forGui.gridR + 4
		allExts = parseFileExtensions(argD.extras)[0]
		theBtn = None
		if len(allExts) > 0:
			theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateMultiFileName(theBox, tkinter.filedialog.asksaveasfilename(title = argD.name, filetypes = tuple([(cv, "*" + cv) for cv in allExts]), defaultextension="*.*" )))
		else:
			theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateMultiFileName(theBox, tkinter.filedialog.asksaveasfilename(title = argD.name )))
		theBtn.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theBtn, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curRetOpts = forData[2].get('1.0', tkinter.END).split("\n")
		curSigOn = 0
		for cv in curRetOpts:
			cvs = cv.strip()
			if len(cvs) == 0:
				continue
			if curSigOn == 0:
				toRet.append(forData[0].sigils[0])
				curSigOn = 1
			if (curSigOn > 1) and (len(cvs) >= 2) and (cvs[0] == '-') and (cvs[1] == '-'):
				toRet.append("--.-")
				toRet.append(forData[0].sigils[0])
				toRet.append(cvs)
				curSigOn = 0
			else:
				toRet.append(cvs)
				curSigOn = 2
		return toRet


class ArgumentFlavorStringFolderRead(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "\\fB" + manSan(argD.sigils[0]) + "\\fR \\fIDIR\\fR"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theBox = tkinter.Entry(forGui.myCanvas)
		theBox.grid(column = 1, row = forGui.gridR)
		theBox.delete(0, tkinter.END)
		initLen = struct.unpack(">Q", argD.extras[0:8])[0]
		initValue = str(argD.extras[8:(initLen+8)], "utf-8")
		theBox.insert(0, initValue)
		theBtn = None
		theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateSingleFileName(theBox, tkinter.filedialog.askdirectory(title = argD.name, mustexist=True)))
		theBtn.grid(column = 2, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theBtn, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curVal = forData[2].get().strip()
		if len(curVal) > 0:
			toRet.append(forData[0].sigils[0])
			toRet.append(curVal)
		return toRet


class ArgumentFlavorStringFolderWrite(ArgumentFlavor):
	def manpageSynopMang(self,argD):
		return "\\fB" + manSan(argD.sigils[0]) + "\\fR \\fIDIR\\fR"
	def makeArgGUI(self,forGui,argD):
		theNam = tkinter.Label(forGui.myCanvas, text = argD.name)
		theNam.grid(column = 0, row = forGui.gridR)
		theBox = tkinter.Entry(forGui.myCanvas)
		theBox.grid(column = 1, row = forGui.gridR)
		theBox.delete(0, tkinter.END)
		initLen = struct.unpack(">Q", argD.extras[0:8])[0]
		initValue = str(argD.extras[8:(initLen+8)], "utf-8")
		theBox.insert(0, initValue)
		theBtn = None
		theBtn = tkinter.Button(forGui.myCanvas, text='Browse', command = lambda: forGui.updateSingleFileName(theBox, tkinter.filedialog.askdirectory(title = argD.name, mustexist=False)))
		theBtn.grid(column = 2, row = forGui.gridR)
		theLab = tkinter.Label(forGui.myCanvas, text = argD.summary)
		theLab.grid(column = 3, row = forGui.gridR)
		forGui.gridR = forGui.gridR + 1
		forGui.argPassFlavs.append(self)
		forGui.argPassCaps.append([argD, theNam, theBox, theBtn, theLab])
	def getGUIArgTexts(self,forGui,forData):
		toRet = []
		curVal = forData[2].get().strip()
		if len(curVal) > 0:
			toRet.append(forData[0].sigils[0])
			toRet.append(curVal)
		return toRet


argFlavorMap = {}
argFlavorMap[("flag","")] = ArgumentFlavorFlag()
argFlavorMap[("enum","")] = ArgumentFlavorEnum()
argFlavorMap[("int","")] = ArgumentFlavorInt()
argFlavorMap[("intvec","")] = ArgumentFlavorIntVec()
argFlavorMap[("intvecg","")] = ArgumentFlavorIntGreedVec()
argFlavorMap[("float","")] = ArgumentFlavorFloat()
argFlavorMap[("floatvec","")] = ArgumentFlavorFloatVec()
argFlavorMap[("floatvecg","")] = ArgumentFlavorFloatGreedVec()
argFlavorMap[("string","")] = ArgumentFlavorString()
argFlavorMap[("stringvec","")] = ArgumentFlavorStringVec()
argFlavorMap[("stringvecg","")] = ArgumentFlavorStringGreedVec()
argFlavorMap[("string","fileread")] = ArgumentFlavorStringFileRead()
argFlavorMap[("stringvec","fileread")] = ArgumentFlavorStringVecFileRead()
argFlavorMap[("stringvecg","fileread")] = ArgumentFlavorStringGreedVecFileRead()
argFlavorMap[("string","filewrite")] = ArgumentFlavorStringFileWrite()
argFlavorMap[("stringvec","filewrite")] = ArgumentFlavorStringVecFileWrite()
argFlavorMap[("stringvecg","filewrite")] = ArgumentFlavorStringGreedVecFileWrite()
argFlavorMap[("string","folderread")] = ArgumentFlavorStringFolderRead()
argFlavorMap[("string","folderwrite")] = ArgumentFlavorStringFolderWrite()


class ArgumentGUIFrame(tkinter.Frame):
	'''Let the user pick arguments for a program.'''
	def updateSingleFileName(self, forBox, newName):
		if newName:
			forBox.delete(0, tkinter.END)
			forBox.insert(0, newName)
	def updateMultiFileName(self, forBox, newName):
		if isinstance(newName, str):
			forBox.insert(tkinter.END, '\n' + newName)
		elif isinstance(newName, tuple):
			for cname in newName:
				forBox.insert(tkinter.END, '\n' + cname)
	def __init__(self,forProg,checkPref,master=None):
		super().__init__(master)
		self.master = master
		self.finalArgs = None
		'''The final arguments: only set on OK.'''
		self.progCheckPref = checkPref
		'''The way to run the program: used to check arguments.'''
		self.mainProg = forProg
		'''The program this is for.'''
		self.argPassFlavs = []
		'''The flavors of the arguments this works with'''
		self.argPassCaps = []
		'''Extra data for the arguments this can work with.'''
		self.myCanvas = tkinter.Canvas(self)
		'''The main canvas of this thing'''
		self.gridR = 0
		'''The grid row this is on.'''
		self.valueCache = {}
		'''Named things to refer back to.'''
		# add a summary line
		self.sumLine = tkinter.Label(self.myCanvas, text = forProg.summary)
		self.sumLine.grid(column = 0, row = self.gridR, columnspan=3)
		self.gridR = self.gridR + 1
		# make stuff for the arguments
		for i in range(len(forProg.arguments)):
			carg = forProg.arguments[i]
			if not carg.isPublic:
				continue
			argTp = (carg.mainFlavor, carg.subFlavor)
			if not (argTp in argFlavorMap):
				continue
			argFlav = argFlavorMap[argTp]
			argFlav.makeArgGUI(self,carg)
		# make the ok and cancel buttons
		self.canBtn = tkinter.Button(self.myCanvas, text = 'Cancel', command = lambda: self.actionCancel())
		self.canBtn.grid(column = 0, row = self.gridR)
		self.goBtn = tkinter.Button(self.myCanvas, text = 'OK', command = lambda: self.actionGo())
		self.goBtn.grid(column = 1, row = self.gridR)
		# wrap it up
		self.myCanvas.pack(side=tkinter.LEFT)
		self.pack()
	def actionCancel(self):
		self.master.destroy()
	def actionGo(self):
		retArgs = []
		# get the arguments
		for i in range(len(self.argPassFlavs)):
			curGArgs = self.argPassFlavs[i].getGUIArgTexts(self,self.argPassCaps[i])
			retArgs.extend(curGArgs)
		# see if the program likes them
		curRunP = subprocess.Popen(self.progCheckPref + ["--help_id10t"] + retArgs,stderr=subprocess.PIPE,stdout=subprocess.DEVNULL,stdin=subprocess.DEVNULL)
		allErrT = curRunP.stderr.read()
		curRunP.stderr.close()
		retCode = curRunP.wait()
		if (len(allErrT) > 0) or (retCode != 0):
			tkinter.messagebox.showerror(title="Problem with arguments!", message=str(allErrT,"utf-8"))
			return
		# all good, set and quit
		self.finalArgs = retArgs
		self.master.destroy()


class ProgramSelectGUIFrame(tkinter.Frame):
	'''Pick a program to run.'''
	def __init__(self,forSet,master=None):
		super().__init__(master)
		self.master = master
		self.finalPI = None
		'''The final selection: only set on OK.'''
		self.myCanvas = tkinter.Canvas(self)
		'''The main canvas of this thing'''
		# make the selections
		gridR = 0
		self.progButs = []
		self.progTexts = []
		for i in range(len(forSet.progNames)):
			progBut = tkinter.Button(self.myCanvas, text = forSet.progNames[i], command = self.helpMakeGo(i))
			progBut.grid(column = 0, row = gridR)
			progTxt = tkinter.Label(self.myCanvas, text = forSet.progSummaries[i])
			progTxt.grid(column = 1, row = gridR)
			self.progButs.append(progBut)
			self.progTexts.append(progTxt)
			gridR = gridR + 1
		# make the cancel button
		self.canBtn = tkinter.Button(self.myCanvas, text = 'Cancel', command = lambda: self.actionCancel())
		self.canBtn.grid(column = 0, row = gridR)
		# make a scroll bar
		self.myCanvas.pack(side=tkinter.LEFT)
		mainscr = tkinter.Scrollbar(self, command=self.myCanvas.yview)
		mainscr.pack(side=tkinter.LEFT, fill='y')
		self.myCanvas.configure(yscrollcommand = mainscr.set)
		self.pack()
	def actionCancel(self):
		self.master.destroy()
	def actionGo(self,forProgI):
		self.finalPI = forProgI
		self.master.destroy()
	def helpMakeGo(self,forInd):
		def doIt():
			self.actionGo(forInd)
		return doIt


def htmlSan(theStr):
	'''
	Make a string safe for html.
	@param theStr: The string to make safe.
	@return: The safe string.
	'''
	outStr = theStr
	outStr = outStr.replace("&","&amp;")
	outStr = outStr.replace('"',"&quot;")
	outStr = outStr.replace("'","&apos;")
	outStr = outStr.replace("<","&lt;")
	outStr = outStr.replace(">","&gt;")
	return outStr


def writeHTMLForProg(progData,toStr):
	'''
	Write an HTML section for a program.
	@param progData: The parsed data for the program.
	@param toStr: The stream to write to.
	'''
	toStr.write("<P>" + htmlSan(progData.usage) + "</P>\n")
	if len(progData.description) > 0:
		allDesc = [cv.strip() for cv in progData.description.split("\n")]
		for cdesc in allDesc:
			toStr.write("<P>" + htmlSan(cdesc) + "</P>\n")
	else:
		toStr.write("<P>" + htmlSan(progData.summary) + "</P>\n")
	toStr.write("<TABLE>\n")
	for carg in progData.arguments:
		if not carg.isPublic:
			continue
		toStr.write("	<TR>\n")
		toStr.write("		<TD>" + htmlSan(carg.name) + "</TD>\n")
		toStr.write("		<TD>" + "/".join([htmlSan(cv) for cv in carg.sigils]) + "</TD>\n")
		toStr.write("		<TD>" + htmlSan(carg.summary) + "</TD>\n")
		toStr.write("	</TR>\n")
		if len(carg.description) > 0:
			toStr.write("	<TR>\n")
			toStr.write("		<TD></TD>\n")
			toStr.write("		<TD></TD>\n")
			allDesc = [cv.strip() for cv in carg.description.split("\n")]
			toStr.write("		<TD>\n")
			for cdesc in allDesc:
				toStr.write("			<P>" + htmlSan(cdesc) + "</P>\n")
			toStr.write("		</TD>\n")
			toStr.write("	</TR>\n")
		if len(carg.usage) > 0:
			toStr.write("	<TR>\n")
			toStr.write("		<TD></TD>\n")
			toStr.write("		<TD></TD>\n")
			toStr.write("		<TD>" + htmlSan(carg.usage) + "</TD>\n")
			toStr.write("	</TR>\n")
	toStr.write("</TABLE>\n")


def writeHTMLProgram(progPath,toStr):
	'''
	Write HTML documentation for a program.
	@param progPath: The path to the program, and any forced arguments.
	@param toStr: The stream to write to.
	'''
	progD = getProgramInfos(progPath,None)
	toStr.write("<!DOCTYPE html>\n")
	toStr.write("<HTML>\n")
	toStr.write("<HEAD>\n")
	toStr.write("	<TITLE>" + htmlSan(progD.name) + "</TITLE>\n")
	toStr.write("</HEAD>\n")
	toStr.write("<BODY>\n")
	toStr.write("<P>" + htmlSan(progD.name) + "</P>\n")
	writeHTMLForProg(progD,toStr)
	toStr.write("</BODY>\n")
	toStr.write("</HTML>\n")


def writeHTMLProgramSet(progPath,toStr):
	'''
	Write HTML documentation for a program set.
	@param progPath: The path to the program, and any forced arguments.
	@param toStr: The stream to write to.
	'''
	progS = getProgramSetInfos(progPath)
	toStr.write("<!DOCTYPE html>\n")
	toStr.write("<HTML>\n")
	toStr.write("<HEAD>\n")
	toStr.write("	<TITLE>" + htmlSan(progS.name) + "</TITLE>\n")
	toStr.write("</HEAD>\n")
	toStr.write("<BODY>\n")
	toStr.write("<P>" + htmlSan(progS.name) + "</P>\n")
	toStr.write("<P>" + htmlSan(progS.summary) + "</P>\n")
	toStr.write("<TABLE>\n")
	allProgDs = []
	for i in range(len(progS.progNames)):
		curNam = progS.progNames[i]
		curSum = progS.progSummaries[i]
		curProD = getProgramInfos(progPath,curNam)
		allProgDs.append(curProD)
		toStr.write("	<TR>\n")
		toStr.write("		<TD><A HREF=\"#prog" + repr(i) + "\">" + htmlSan(curNam) + "</A></TD>\n")
		toStr.write("		<TD>" + htmlSan(curSum) + "</TD>\n")
		toStr.write("	</TR>\n")
	toStr.write("</TABLE>\n")
	for i in range(len(progS.progNames)):
		toStr.write("<HR>\n")
		toStr.write("<P><A ID=\"prog" + repr(i) + "\">" + htmlSan(progS.name) + " " + htmlSan(progS.progNames[i]) + "</A></P>\n")
		writeHTMLForProg(allProgDs[i],toStr)
	toStr.write("</BODY>\n")
	toStr.write("</HTML>\n")


def manSan(theStr):
	'''
	Sanitize a string for a manpage.
	'''
	outStr = theStr
	outStr = outStr.replace("-", "\\-")
	return outStr


def writeManpageForProg(psetName,progData,toStr):
	'''
	Write a manpage for a program.
	@param psetName: The name of the program set it is in, if any.
	@param progData: The parsed data for the program.
	@param toStr: The stream to write to.
	'''
	workPName = progData.name
	if len(psetName) > 0:
		workPName = psetName + " " + workPName
		toStr.write(".TH " + manSan(psetName).upper() + "-" + manSan(progData.name).upper() + " 1\n")
	else:
		toStr.write(".TH " + manSan(progData.name).upper() + " 1\n")
	toStr.write(".SH NAME\n")
	toStr.write(manSan(workPName) + " \\- " + manSan(progData.summary) + "\n")
	toStr.write(".SH SYNOPSIS\n")
	toStr.write(".B " + manSan(workPName) + "\n")
	for carg in progData.arguments:
		if not carg.isPublic:
			continue
		argTp = (carg.mainFlavor, carg.subFlavor)
		if not (argTp in argFlavorMap):
			continue
		argFlav = argFlavorMap[argTp]
		toStr.write(argFlav.manpageSynopMang(carg) + "\n")
	toStr.write(".SH DESCRIPTION\n")
	toStr.write(".B " + manSan(workPName) + "\n")
	if len(progData.description) > 0:
		allDesc = [cv.strip() for cv in progData.description.split("\n")]
		for cdesc in allDesc:
			toStr.write(manSan(cdesc) + "\n")
	else:
		toStr.write(manSan(progData.summary) + "\n")
	toStr.write(".SH OPTIONS\n")
	for carg in progData.arguments:
		argTp = (carg.mainFlavor, carg.subFlavor)
		if not (argTp in argFlavorMap):
			continue
		argFlav = argFlavorMap[argTp]
		toStr.write(".TP\n")
		toStr.write(argFlav.manpageSynopMang(carg) + "\n")
		toStr.write(manSan(carg.summary) + "\n")


def writeManpageProgram(progPath,toStr):
	'''
	Write manpage documentation for a program.
	@param progPath: The path to the program.
	@param subPName: The sub-name of the program of interest.
	@param toStr: The stream to write to.
	'''
	progD = getProgramInfos(progPath,None)
	writeManpageForProg("",progD,toStr)


def writeManpageProgramSet(progPath,toFold,prefix):
	'''
	Write HTML documentation for a program set.
	@param progPath: The path to the program.
	@param toFold: The folder to write to.
	@param prefix: The prefix for the file names.
	'''
	progS = getProgramSetInfos(progPath)
	toStr = open(os.path.join(toFold,prefix + ".1"), "w")
	toStr.write(".TH " + manSan(progS.name).upper() + " 1\n")
	toStr.write(".SH NAME\n")
	toStr.write(manSan(progS.name) + " \\- " + manSan(progS.summary) + "\n")
	toStr.write(".SH SYNOPSIS\n")
	toStr.write(".B " + manSan(progS.name) + "\n")
	toStr.write("\\fIPROGRAM\\fR \\fIARGS\\fR ...\n")
	toStr.write(".SH DESCRIPTION\n")
	toStr.write(".B " + manSan(progS.name) + "\n")
	toStr.write(manSan(progS.summary) + "\n")
	toStr.write(".SH OPTIONS\n")
	allProgDs = []
	for i in range(len(progS.progNames)):
		curNam = progS.progNames[i]
		curSum = progS.progSummaries[i]
		curProD = getProgramInfos(progPath,curNam)
		allProgDs.append(curProD)
		toStr.write(".TP\n")
		toStr.write("\\fB" + manSan(progS.name) + "\\fR " + "\\fB" + manSan(curNam) + "\\fR \\fIARGS\\fR ...\n")
		toStr.write(manSan(curSum) + "\n")
	toStr.close()
	for i in range(len(allProgDs)):
		subStr = open(os.path.join(toFold,prefix + "-" + progS.progNames[i] + ".1"), "w")
		writeManpageForProg(progS.name,allProgDs[i],subStr)
		subStr.close()


def mdSan(toSan):
	'''
	Sanitize some text for markdown.
	@param toSan: The text to sanitize.
	@return: The sanitized text.
	'''
	#\`*_{}[]()#+-.!
	outStr = toSan
	outStr = outStr.replace("\\", "\\\\")
	outStr = outStr.replace("`", "\\`")
	outStr = outStr.replace("*", "\\*")
	outStr = outStr.replace("_", "\\_")
	outStr = outStr.replace("{", "\\{")
	outStr = outStr.replace("}", "\\}")
	outStr = outStr.replace("[", "\\[")
	outStr = outStr.replace("]", "\\]")
	outStr = outStr.replace("(", "\\(")
	outStr = outStr.replace(")", "\\)")
	outStr = outStr.replace("#", "\\#")
	outStr = outStr.replace("+", "\\+")
	outStr = outStr.replace("-", "\\-")
	outStr = outStr.replace(".", "\\.")
	outStr = outStr.replace("!", "\\!")
	return outStr


def writeMarkdownForProg(progData,toStr):
	'''
	Write an markdown section for a program.
	@param progData: The parsed data for the program.
	@param toStr: The stream to write to.
	'''
	toStr.write("# " + mdSan(progData.name) + "\n\n")
	if len(progData.description) > 0:
		allDesc = [cv.strip() for cv in progData.description.split("\n")]
		for cdesc in allDesc:
			toStr.write(mdSan(cdesc.strip()) + "\n\n")
	else:
		toStr.write(mdSan(progData.summary.strip()) + "\n\n")
	for carg in progData.arguments:
		if not carg.isPublic:
			continue
		toStr.write("- ")
		if len(carg.usage.strip()) > 0:
			toStr.write(mdSan(carg.usage.strip()))
		elif len(carg.sigils) > 0:
			toStr.write(mdSan(carg.sigils[0].strip()))
		else:
			toStr.write(mdSan(carg.name.strip()))
		toStr.write("\n\n")
		if len(carg.description) > 0:
			allDesc = [cv.strip() for cv in carg.description.split("\n")]
			for cdesc in allDesc:
				toStr.write("    " + mdSan(cdesc.strip()) + "\n")
		else:
			toStr.write("    " + mdSan(carg.summary.strip()) + "\n")


def writeMarkdownProgram(progPath,toStr):
	'''
	Write markdown documentation for a program.
	@param progPath: The path to the program.
	@param toStr: The stream to write to.
	'''
	progD = getProgramInfos(progPath,None)
	writeMarkdownForProg(progD,toStr)


def writeMarkdownProgramSet(progPath,toStr):
	'''
	Write markdown documentation for a program set.
	@param progPath: The path to the program.
	@param toStr: The stream to write to.
	'''
	progS = getProgramSetInfos(progPath)
	toStr.write("# " + mdSan(progS.name.strip()) + "\n\n")
	toStr.write(mdSan(progS.summary.strip()) + "\n\n")
	allProgDs = []
	for i in range(len(progS.progNames)):
		curNam = progS.progNames[i]
		curSum = progS.progSummaries[i]
		curProD = getProgramInfos(progPath,curNam)
		allProgDs.append(curProD)
		toStr.write("- [" + mdSan(curNam.strip()) + "](#" + curNam.strip().replace(" ","-").lower() + ") : " + mdSan(curSum.strip()) + "\n")
	toStr.write("\n\n")
	for i in range(len(progS.progNames)):
		writeMarkdownForProg(allProgDs[i],toStr)
		toStr.write("\n\n")


def latexSan(toSan):
	'''
	Sanitize a string for latex consumption.
	@param toSan: The string to sanitize.
	@return: The sanitized string.
	'''
	outStr = toSan
	outStr = outStr.replace("\\", "\\textbackslash ")
	outStr = outStr.replace("&", "\\&")
	outStr = outStr.replace("%", "\\%")
	outStr = outStr.replace("$", "\\$")
	outStr = outStr.replace("#", "\\#")
	outStr = outStr.replace("_", "\\_")
	outStr = outStr.replace("{", "\\{")
	outStr = outStr.replace("}", "\\}")
	outStr = outStr.replace("~", "\\textasciitilde ")
	outStr = outStr.replace("^", "\\textasciicircum ")
	outStr = outStr.replace("--", "-{}-")
	return outStr


def writeLatexForProg(progData,toStr):
	'''
	Write an latex section for a program.
	@param progData: The parsed data for the program.
	@param toStr: The stream to write to.
	'''
	toStr.write("\\section{" + latexSan(progData.name) + "}\n")
	if len(progData.description) > 0:
		allDesc = [cv.strip() for cv in progData.description.split("\n")]
		for cdesc in allDesc:
			toStr.write(latexSan(cdesc.strip()) + "\n\n")
	else:
		toStr.write(latexSan(progData.summary.strip()) + "\n\n")
	for carg in progData.arguments:
		if not carg.isPublic:
			continue
		toStr.write("\\textbf{")
		if len(carg.sigils) > 0:
			toStr.write(latexSan(carg.sigils[0].strip()))
		else:
			toStr.write(latexSan(carg.name.strip()))
		toStr.write("}\n\n")
		if len(carg.usage) > 0:
			toStr.write(latexSan(carg.usage.strip()) + "\n\n")
		if len(carg.description) > 0:
			allDesc = [cv.strip() for cv in carg.description.split("\n")]
			for cdesc in allDesc:
				toStr.write(latexSan(cdesc.strip()) + "\n")
		else:
			toStr.write(latexSan(carg.summary.strip()) + "\n")
		toStr.write("\n\n")


def writeLatexProgram(progPath,toStr):
	'''
	Write latex documentation for a program.
	@param progPath: The path to the program.
	@param toStr: The stream to write to.
	'''
	progD = getProgramInfos(progPath,None)
	toStr.write("\\documentclass{article}\n")
	toStr.write("\\title{" + latexSan(progD.name.strip()) + "}\n")
	toStr.write("\\begin{document}\n")
	writeLatexForProg(progD,toStr)
	toStr.write("\\end{document}\n")


def writeLatexProgramSet(progPath,toStr):
	'''
	Write latex documentation for a program set.
	@param progPath: The path to the program.
	@param toStr: The stream to write to.
	'''
	progS = getProgramSetInfos(progPath)
	toStr.write("\\documentclass{article}\n")
	toStr.write("\\title{" + latexSan(progS.name.strip()) + "}\n")
	toStr.write("\\begin{document}\n")
	toStr.write(latexSan(progS.name.strip()) + "\n\n")
	toStr.write(latexSan(progS.summary.strip()) + "\n\n")
	toStr.write("\\tableofcontents\n\n")
	for i in range(len(progS.progNames)):
		curNam = progS.progNames[i]
		curSum = progS.progSummaries[i]
		curProD = getProgramInfos(progPath,curNam)
		writeLatexForProg(curProD,toStr)
		toStr.write("\n\n")
	toStr.write("\\end{document}\n")


def runGuiProgram(progPath,subPName = None):
	'''
	Run a gui for a specific program.
	@param progPath: The path to the program.
	@param subPName: The name of the sub-program to run, if relevant.
	'''
	progD = getProgramInfos(progPath,subPName)
	# get the arguments
	root = tkinter.Tk()
	root.title(progD.name)
	checkPref = progPath[:]
	if not (subPName is None):
		checkPref.append(subPName)
	app = ArgumentGUIFrame(progD, checkPref, root)
	app.mainloop()
	if app.finalArgs is None:
		return
	fullCallAs = progPath[:]
	if not (subPName is None):
		fullCallAs.append(subPName)
	fullCallAs.extend(app.finalArgs)
	print(" ".join(fullCallAs))
	# run the thing
	curRunP = subprocess.Popen(fullCallAs,stdout=subprocess.DEVNULL,stdin=subprocess.DEVNULL,stderr=subprocess.PIPE)
	allErrT = curRunP.stderr.read()
	allErrT = str(allErrT, "utf-8").strip()
	curRunP.stderr.close()
	if (curRunP.wait() != 0) or (len(allErrT) != 0):
		tkinter.messagebox.showerror(title="Problem running program!", message=allErrT)
		raise IOError("Problem running program.")


def runGuiProgramSet(progPath):
	'''
	Run a gui for a program set.
	@param progPath: The path to the program.
	'''
	progS = getProgramSetInfos(progPath)
	# select the thing to run
	root = tkinter.Tk()
	root.title(progS.name)
	app = ProgramSelectGUIFrame(progS, root)
	app.mainloop()
	if app.finalPI is None:
		return
	# select the arguments for the thing
	runGuiProgram(progPath, progS.progNames[app.finalPI])


def commonProgSetup(forPro):
	forPro.version = "ArgMang 0.0\nCopyright (C) 2023 Benjamin Crysup\nLicense LGPLv3: GNU LGPL version 3\nThis is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n"
	# a string option
	forPro.strOpt = whodunargs.ArgumentOptionString("--prog","Executable","The executable to run to get data.","--prog python3")
	forPro.options.append(forPro.strOpt)
	# any arguments to the thing
	forPro.strvOpt = whodunargs.ArgumentOptionStringVector("--arg","Required Arguments","Any required arguments to the program.","--arg WDAExample.py")
	forPro.options.append(forPro.strvOpt)


class DocumentHTMLProgram(whodunargs.StandardProgram):
	'''Document to an html file.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "htmlp"
		self.summary = "Generate html documentation for a program."
		self.usage = "python3 htmlp --prog EXE --arg ARG --out FILE"
		commonProgSetup(self)
		# where to write the documentation: if not specified, stdout
		self.filewOpt = whodunargs.ArgumentOptionFileWrite("--out","Output","Where to write the html.","--out /home/ben/WDAExample.html")
		self.filewOpt.validExts.add(".html")
		self.options.append(self.filewOpt)
	def idiotCheck(self):
		if len(self.strOpt.value) == 0:
			raise ValueError("No program provided.")
	def baseRun(self):
		toOut = self.useOut
		if len(self.filewOpt.value) > 0:
			toOut = open(self.filewOpt.value,"wb")
		toOut = io.TextIOWrapper(toOut, encoding="utf-8")
		writeHTMLProgram([self.strOpt.value] + self.strvOpt.value,toOut)
		toOut.close()


class DocumentSetHTMLProgram(whodunargs.StandardProgram):
	'''Document to an html file.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "htmls"
		self.summary = "Generate html documentation for a program set."
		self.usage = "python3 htmls --prog EXE --arg ARG --out FILE"
		commonProgSetup(self)
		# where to write the documentation: if not specified, stdout
		self.filewOpt = whodunargs.ArgumentOptionFileWrite("--out","Output","Where to write the html.","--out /home/ben/ArgMang.html")
		self.filewOpt.validExts.add(".html")
		self.options.append(self.filewOpt)
	def idiotCheck(self):
		if len(self.strOpt.value) == 0:
			raise ValueError("No program provided.")
	def baseRun(self):
		toOut = self.useOut
		if len(self.filewOpt.value) > 0:
			toOut = open(self.filewOpt.value,"wb")
		toOut = io.TextIOWrapper(toOut, encoding="utf-8")
		writeHTMLProgramSet([self.strOpt.value] + self.strvOpt.value,toOut)
		toOut.close()


class DocumentManpageProgram(whodunargs.StandardProgram):
	'''Document to a troff file.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "manp"
		self.summary = "Generate a manpage for a program."
		self.usage = "python3 manp --prog EXE --arg ARG --out FILE"
		commonProgSetup(self)
		# where to write the documentation: if not specified, stdout
		self.filewOpt = whodunargs.ArgumentOptionFileWrite("--out","Output","Where to write the troff file.","--out /home/ben/WDAExample.1")
		self.options.append(self.filewOpt)
	def idiotCheck(self):
		if len(self.strOpt.value) == 0:
			raise ValueError("No program provided.")
	def baseRun(self):
		toOut = self.useOut
		if len(self.filewOpt.value) > 0:
			toOut = open(self.filewOpt.value,"wb")
		toOut = io.TextIOWrapper(toOut, encoding="utf-8")
		writeManpageProgram([self.strOpt.value] + self.strvOpt.value,toOut)
		toOut.close()


class DocumentSetManpageProgram(whodunargs.StandardProgram):
	'''Document to troff files.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "mans"
		self.summary = "Generate manpages for a program set."
		self.usage = "python3 mans --prog EXE --arg ARG --out FOLDER"
		commonProgSetup(self)
		# where to write the documentation
		self.filewOpt = whodunargs.ArgumentOptionFolderWrite("--out","Output","Where to write the files.","--out /home/ben/ArgMang/")
		self.options.append(self.filewOpt)
		# the prefix for the output name
		self.prefOpt = whodunargs.ArgumentOptionString("--pref","Prefix","The prefix for the file names.","--pref argmang")
		self.options.append(self.prefOpt)
	def idiotCheck(self):
		if len(self.strOpt.value) == 0:
			raise ValueError("No program provided.")
		if len(self.prefOpt.value) == 0:
			raise ValueError("No prefix provided.")
	def baseRun(self):
		toOut = "."
		if len(self.filewOpt.value) > 0:
			toOut = self.filewOpt.value
		writeManpageProgramSet([self.strOpt.value] + self.strvOpt.value,toOut,self.prefOpt.value)


class DocumentMarkdownProgram(whodunargs.StandardProgram):
	'''Document to an md file.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "mdp"
		self.summary = "Generate markdown documentation for a program."
		self.usage = "python3 mdp --prog EXE --arg ARG --out FILE"
		commonProgSetup(self)
		# where to write the documentation: if not specified, stdout
		self.filewOpt = whodunargs.ArgumentOptionFileWrite("--out","Output","Where to write the md file.","--out /home/ben/WDAExample.md")
		self.filewOpt.validExts.add(".md")
		self.options.append(self.filewOpt)
	def idiotCheck(self):
		if len(self.strOpt.value) == 0:
			raise ValueError("No program provided.")
	def baseRun(self):
		toOut = self.useOut
		if len(self.filewOpt.value) > 0:
			toOut = open(self.filewOpt.value,"wb")
		toOut = io.TextIOWrapper(toOut, encoding="utf-8")
		writeMarkdownProgram([self.strOpt.value] + self.strvOpt.value,toOut)
		toOut.close()


class DocumentSetMarkdownProgram(whodunargs.StandardProgram):
	'''Document to an md file.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "mds"
		self.summary = "Generate markdown documentation for a program set."
		self.usage = "python3 mds --prog EXE --arg ARG --out FILE"
		commonProgSetup(self)
		# where to write the documentation: if not specified, stdout
		self.filewOpt = whodunargs.ArgumentOptionFileWrite("--out","Output","Where to write the md file.","--out /home/ben/ArgMang.md")
		self.filewOpt.validExts.add(".md")
		self.options.append(self.filewOpt)
	def idiotCheck(self):
		if len(self.strOpt.value) == 0:
			raise ValueError("No program provided.")
	def baseRun(self):
		toOut = self.useOut
		if len(self.filewOpt.value) > 0:
			toOut = open(self.filewOpt.value,"wb")
		toOut = io.TextIOWrapper(toOut, encoding="utf-8")
		writeMarkdownProgramSet([self.strOpt.value] + self.strvOpt.value,toOut)
		toOut.close()


class DocumentLatexProgram(whodunargs.StandardProgram):
	'''Document to an latex file.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "latexp"
		self.summary = "Generate latex documentation for a program."
		self.usage = "python3 latexp --prog EXE --arg ARG --out FILE"
		commonProgSetup(self)
		# where to write the documentation: if not specified, stdout
		self.filewOpt = whodunargs.ArgumentOptionFileWrite("--out","Output","Where to write the latex file.","--out /home/ben/WDAExample.tex")
		self.filewOpt.validExts.add(".tex")
		self.options.append(self.filewOpt)
	def idiotCheck(self):
		if len(self.strOpt.value) == 0:
			raise ValueError("No program provided.")
	def baseRun(self):
		toOut = self.useOut
		if len(self.filewOpt.value) > 0:
			toOut = open(self.filewOpt.value,"wb")
		toOut = io.TextIOWrapper(toOut, encoding="utf-8")
		writeLatexProgram([self.strOpt.value] + self.strvOpt.value,toOut)
		toOut.close()


class DocumentSetLatexProgram(whodunargs.StandardProgram):
	'''Document to a latex file.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "latexs"
		self.summary = "Generate latex documentation for a program set."
		self.usage = "python3 latexs --prog EXE --arg ARG --out FILE"
		commonProgSetup(self)
		# where to write the documentation: if not specified, stdout
		self.filewOpt = whodunargs.ArgumentOptionFileWrite("--out","Output","Where to write the latex file.","--out /home/ben/ArgMang.tex")
		self.filewOpt.validExts.add(".tex")
		self.options.append(self.filewOpt)
	def idiotCheck(self):
		if len(self.strOpt.value) == 0:
			raise ValueError("No program provided.")
	def baseRun(self):
		toOut = self.useOut
		if len(self.filewOpt.value) > 0:
			toOut = open(self.filewOpt.value,"wb")
		toOut = io.TextIOWrapper(toOut, encoding="utf-8")
		writeLatexProgramSet([self.strOpt.value] + self.strvOpt.value,toOut)
		toOut.close()


class GuiProgram(whodunargs.StandardProgram):
	'''Run a gui for a program.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "guip"
		self.summary = "Show a helpful gui for a program."
		self.usage = "python3 guip --prog EXE --arg ARG"
		self.description = "Build a gui to pick arguments for a program, then run.\nMake sure your version of python has a working version of tkinter."
		commonProgSetup(self)
	def idiotCheck(self):
		if len(self.strOpt.value) == 0:
			raise ValueError("No program provided.")
	def baseRun(self):
		runGuiProgram([self.strOpt.value] + self.strvOpt.value,None)


class GuiSetProgram(whodunargs.StandardProgram):
	'''Run a gui for a program set.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "guis"
		self.summary = "Show a helpful gui for a program set."
		self.description = "Build a gui to pick arguments for a set of programs, then run.\nMake sure your version of python has a working version of tkinter."
		self.usage = "python3 guis --prog EXE --arg ARG"
		commonProgSetup(self)
	def idiotCheck(self):
		if len(self.strOpt.value) == 0:
			raise ValueError("No program provided.")
	def baseRun(self):
		runGuiProgramSet([self.strOpt.value] + self.strvOpt.value)


class ArgMangProgramSet(whodunargs.StandardProgramSet):
	def __init__(self):
		whodunargs.StandardProgramSet.__init__(self)
		self.name = "ArgMang"
		self.summary = "Do things to help with arguments."
		self.version = "ArgMang 0.0\nCopyright (C) 2023 Benjamin Crysup\nLicense LGPLv3: GNU LGPL version 3\nThis is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n"
		self.programs["htmlp"] = lambda : DocumentHTMLProgram()
		self.programs["htmls"] = lambda : DocumentSetHTMLProgram()
		self.programs["manp"] = lambda : DocumentManpageProgram()
		self.programs["mans"] = lambda : DocumentSetManpageProgram()
		self.programs["mdp"] = lambda : DocumentMarkdownProgram()
		self.programs["mds"] = lambda : DocumentSetMarkdownProgram()
		self.programs["latexp"] = lambda : DocumentLatexProgram()
		self.programs["latexs"] = lambda : DocumentSetLatexProgram()
		self.programs["guip"] = lambda : GuiProgram()
		self.programs["guis"] = lambda : GuiSetProgram()


if __name__ == "__main__":
	toRun = ArgMangProgramSet()
	toRun = toRun.parseArguments(sys.argv[1:])
	if not (toRun is None):
		toRun.run()




