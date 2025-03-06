'''Do things for backups'''

import os
import sys
import time
import shutil
import datetime
import subprocess
import multiprocessing

import tapecom
import whodunargs

commonLicense = "TapeDeck 0.0\nCopyright (C) 2024 Benjamin Crysup\nLicense LGPLv3: GNU LGPL version 3\nThis is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n"


def _openOutput(outFN,baseOut):
	if outFN.strip() != "":
		return open(outFN, "wb")
	return baseOut


def _closeOutput(outFN,outF):
	if outFN.strip() != "":
		outF.close()


class InitDatabaseProgram(whodunargs.StandardProgram):
	'''Initialize an empty database.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "InitDatabase"
		self.summary = "Make a new database."
		self.usage = "python3 TapeDeck InitDatabase"
		self.version = commonLicense
		self.dataRoot = whodunargs.ArgumentOptionFolderWrite("--database","Database Directory","The folder containing the database.","--database C:\\Database\\")
		self.options.append(self.dataRoot)
	def baseRun(self):
		# "load" the database
		tapeDB = tapecom.TapeDatabase(self.dataRoot.value)


class ListFoldersProgram(whodunargs.StandardProgram):
	'''List known folders in a database.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "ListFolders"
		self.summary = "List known data folders."
		self.usage = "python3 TapeDeck ListFolders"
		self.version = commonLicense
		self.dataRoot = whodunargs.ArgumentOptionFolderRead("--database","Database Directory","The folder containing the database.","--database C:\\Database\\")
		self.reportFN = whodunargs.ArgumentOptionFileWrite("--out","Output","The place to write the list.","--out /home/ben/report.txt")
		self.options.append(self.dataRoot)
		self.options.append(self.reportFN)
	def baseRun(self):
		# load the database
		tapeDB = tapecom.TapeDatabase(self.dataRoot.value)
		# open the report
		outF = _openOutput(self.reportFN.value,self.useOut)
		# run down the database
		allFolds = tapeDB.getFolders()
		for cf in allFolds:
			curFD = tapeDB.getFolderData(cf)
			outF.write(bytes(cf,"utf-8"))
			for keyN in curFD:
				outF.write(bytes("\t","utf-8"))
				outF.write(bytes(keyN,"utf-8"))
				outF.write(bytes(":","utf-8"))
				outF.write(bytes(curFD[keyN],"utf-8"))
			outF.write(bytes("\n","utf-8"))
		# and close
		_closeOutput(self.reportFN.value,outF)


class HuntNewProgram(whodunargs.StandardProgram):
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "HuntNewFolders"
		self.summary = "Hunt newly added data folders."
		self.usage = "python3 TapeDeck HuntNewFolders"
		self.version = commonLicense
		self.dataRoot = whodunargs.ArgumentOptionFolderWrite("--database","Database Directory","The folder containing the database.","--database C:\\Database\\")
		self.foldRoot = whodunargs.ArgumentOptionFolderRead("--backfold","Backup Folder","The folder to back up.","--backfold C:\\WorkingSet\\")
		self.reportFN = whodunargs.ArgumentOptionFileWrite("--out","Output","The place to list newly found folders.","--out /home/ben/report.txt")
		self.printOnly = whodunargs.ArgumentOptionFlag("--dry-run","List Only","Just list new folders, do not add to database.")
		self.options.append(self.dataRoot)
		self.options.append(self.foldRoot)
		self.options.append(self.reportFN)
		self.options.append(self.printOnly)
	def baseRun(self):
		# note the current time
		addTime = int(time.time())
		# load the database
		tapeDB = tapecom.TapeDatabase(self.dataRoot.value)
		# and the staging ground
		stageF = tapecom.StagingFolder(self.foldRoot.value)
		# open the output file
		outF = _openOutput(self.reportFN.value,self.useOut)
		# get/add new things
		addedFs = None
		if self.printOnly.value:
			addedFs = stageF.getUnknown(tapeDB)
		else:
			addedFs = stageF.addUnknown(tapeDB)
		# report
		for cf in addedFs:
			outF.write(bytes(cf,"utf-8"))
			outF.write(bytes("\n","utf-8"))
		# close the file
		_closeOutput(self.reportFN.value,outF)


class VerifyFoldersProgram(whodunargs.StandardProgram):
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "VerifyFolders"
		self.summary = "Check that data folders match backup."
		self.usage = "python3 TapeDeck VerifyFolders"
		self.version = commonLicense
		self.dataRoot = whodunargs.ArgumentOptionFolderRead("--database","Database Directory","The folder containing the database.","--database C:\\Database\\")
		self.foldRoot = whodunargs.ArgumentOptionFolderRead("--backfold","Backup Folder","The folder to back up.","--backfold C:\\WorkingSet\\")
		self.reportFN = whodunargs.ArgumentOptionFileWrite("--out","Output","The place to list newly found folders.","--out /home/ben/report.txt")
		self.limitFNs = whodunargs.ArgumentOptionStringGreedyVector("--fold","Limit Folders","The folders to look at. By default, all.","--fold a b")
		self.options.append(self.dataRoot)
		self.options.append(self.foldRoot)
		self.options.append(self.reportFN)
		self.options.append(self.limitFNs)
	def baseRun(self):
		# load the database
		tapeDB = tapecom.TapeDatabase(self.dataRoot.value)
		# and the staging ground
		stageF = tapecom.StagingFolder(self.foldRoot.value)
		# determine which folders to look at
		lookFolds = self.limitFNs.value
		if len(lookFolds) == 0:
			lookFolds = stageF.getKnown(tapeDB)
		# open the output file
		outF = _openOutput(self.reportFN.value,self.useOut)
		# walk through, getting problems
		for cf in lookFolds:
			curPS = stageF.checkFolderHash(tapeDB,cf)
			for cp in curPS:
				outF.write(bytes(cp[0],"utf-8"))
				outF.write(bytes(": ","utf-8"))
				outF.write(bytes(cp[1],"utf-8"))
				outF.write(bytes("\n","utf-8"))
		# close the file
		_closeOutput(self.reportFN.value,outF)


class ListTapesProgram(whodunargs.StandardProgram):
	'''List tapes in a database.'''
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "ListTapes"
		self.summary = "List tapes."
		self.usage = "python3 TapeDeck ListTapes"
		self.version = commonLicense
		self.dataRoot = whodunargs.ArgumentOptionFolderRead("--database","Database Directory","The folder containing the database.","--database C:\\Database\\")
		self.reportFN = whodunargs.ArgumentOptionFileWrite("--out","Output","The place to write the list.","--out /home/ben/report.txt")
		self.options.append(self.dataRoot)
		self.options.append(self.reportFN)
	def baseRun(self):
		# load the database
		tapeDB = tapecom.TapeDatabase(self.dataRoot.value)
		# open the report
		outF = _openOutput(self.reportFN.value,self.useOut)
		# run down the database
		allFolds = tapeDB.getTapes()
		for cf in allFolds:
			curFD = tapeDB.getTapeData(cf)
			outF.write(bytes(cf,"utf-8"))
			for keyN in curFD:
				outF.write(bytes("\t","utf-8"))
				outF.write(bytes(keyN,"utf-8"))
				outF.write(bytes(":","utf-8"))
				outF.write(bytes(curFD[keyN],"utf-8"))
			outF.write(bytes("\n","utf-8"))
		# and close
		_closeOutput(self.reportFN.value,outF)


class ListUntapedProgram(whodunargs.StandardProgram):
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "ListUntaped"
		self.summary = "List folders that are not currently taped."
		self.usage = "python3 TapeDeck ListUntaped"
		self.version = commonLicense
		self.dataRoot = whodunargs.ArgumentOptionFolderWrite("--database","Database Directory","The folder containing the database.","--database C:\\Database\\")
		self.foldRoot = whodunargs.ArgumentOptionFolderRead("--backfold","Backup Folder","The folder to back up.","--backfold C:\\WorkingSet\\")
		self.reportFN = whodunargs.ArgumentOptionFileWrite("--out","Output","The place to list newly found folders.","--out /home/ben/report.txt")
		self.options.append(self.dataRoot)
		self.options.append(self.foldRoot)
		self.options.append(self.reportFN)
	def baseRun(self):
		# load the database
		tapeDB = tapecom.TapeDatabase(self.dataRoot.value)
		# figure out how to get the set of untaped things
		utSet = None
		if len(self.foldRoot.value.strip()) == 0:
			# get all untaped
			allFold = set(tapeDB.getFolders())
			backMap = set(tapeDB.getTapeInverseContentMap())
			utSet = list(allFold.difference(backMap))
		else:
			# get untaped in staging folder
			stageF = tapecom.StagingFolder(self.foldRoot.value)
			utSet = stageF.getUntaped(tapeDB)
		# open the output file
		outF = _openOutput(self.reportFN.value,self.useOut)
		# report
		for cf in utSet:
			curFD = tapeDB.getFolderData(cf)
			outF.write(bytes(cf,"utf-8"))
			for keyN in curFD:
				outF.write(bytes("\t","utf-8"))
				outF.write(bytes(keyN,"utf-8"))
				outF.write(bytes(":","utf-8"))
				outF.write(bytes(curFD[keyN],"utf-8"))
			outF.write(bytes("\n","utf-8"))
		# close the file
		_closeOutput(self.reportFN.value,outF)


class NewTapeProgram(whodunargs.StandardProgram):
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "NewTape"
		self.summary = "Create a new tape."
		self.usage = "python3 TapeDeck NewTape"
		self.version = commonLicense
		self.dataRoot = whodunargs.ArgumentOptionFolderWrite("--database","Database Directory","The folder containing the database.","--database C:\\Database\\")
		self.foldRoot = whodunargs.ArgumentOptionFolderRead("--backfold","Backup Folder","The folder to back up.","--backfold C:\\WorkingSet\\")
		self.tempRoot = whodunargs.ArgumentOptionFolderWrite("--scratch","Scratch Directory","A folder for temporary tar files.","--scratch /tmp/vhs")
		self.reportFN = whodunargs.ArgumentOptionFileWrite("--out","Output","The place to write the final tar file.","--out /dev/st0")
		self.minTotFS = whodunargs.ArgumentOptionInteger("--min","Minimum Size","The minimum number of bytes to push to tape.","--min 9000000000000")
		self.maxTotFS = whodunargs.ArgumentOptionInteger("--max","Maximum Size","The maximum number of bytes to push to tape.","--max 11000000000000")
		self.tapeName = whodunargs.ArgumentOptionString("--name","Tape Name","The name to assign to the tape.","--name T2024_02_07")
		self.options.append(self.dataRoot)
		self.options.append(self.foldRoot)
		self.options.append(self.tempRoot)
		self.options.append(self.reportFN)
		self.options.append(self.minTotFS)
		self.options.append(self.maxTotFS)
		self.options.append(self.tapeName)
		self.minTotFS.value = 0
		self.maxTotFS.value = 11000000000000
	def baseRun(self):
		if len(self.reportFN.value.strip()) == 0:
			raise ValueError("NewTape cannout output data on stdout.")
		# load the database
		tapeDB = tapecom.TapeDatabase(self.dataRoot.value)
		# if tape name already in, complain
		if self.tapeName.value in tapeDB.getTapes():
			raise IOError(self.tapeName.value + " is already a tape.")
		# and the staging ground
		stageF = tapecom.StagingFolder(self.foldRoot.value)
		# note what is not taped
		utSet = stageF.getUntaped(tapeDB)
		utSize = [int(tapeDB.getFolderData(cu)["Size"]) for cu in utSet]
		ntSet = []
		totNtSize = 0
		# if there is not enough to at least make minimum, quit
		if sum(utSize) < self.minTotFS.value:
			raise IOError("Not enough data currently present to justify a new tape.")
		# make the scratch folder
		os.makedirs(self.tempRoot.value, exist_ok=False)
		try:
			# tar up folders until enough is had
			while True:
				oldNTL = len(ntSet)
				totNtSize = self.grabMoreFolds(utSet, utSize, ntSet, totNtSize, tapeDB)
				if len(ntSet) == oldNTL:
					break
			# tar up the pieces into a super-tar
			ntPaths = [cfn + ".tar.gz.aes" for cfn in ntSet]
			spTar = subprocess.Popen(["tar","cf",os.path.abspath(self.reportFN.value)] + ntPaths,stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,cwd=self.tempRoot.value)
			if spTar.wait() != 0:
				raise IOError("Problem writing " + self.reportFN.value)
		except:
			shutil.rmtree(self.tempRoot.value)
			raise
		# delete scratch
		shutil.rmtree(self.tempRoot.value)
		# set up the map
		datMap = {}
		datMap["Size"] = repr(totNtSize)
		datMap["Added"] = repr(int(time.time()))
		# add the tape
		tapeDB.addTape(self.tapeName.value,datMap,ntSet)
	def grabMoreFolds(self,utSet,utSize,ntSet,startSize,tapeDB):
		'''
		Prepare more folders for addition to a tape.
		@param utSet: The list of folders waiting for a tape.
		@param utSize: The size of each folder waiting for a tape.
		@param ntSet: The list of folders ready to put on tape.
		@param startSize: The number of bytes ntSet represents on call.
		@param tapeDB: The tape database.
		@return: The number of bytes ntSet represents on return.
		'''
		maxOneGo = multiprocessing.cpu_count()
		# figure out how many to get in one round (and how many can't even)
		grabInds = []
		putEndSize = startSize
		for pi in range(len(utSet)):
			curUTN = utSet[pi]
			curUTSize = utSize[pi]
			if (putEndSize + curUTSize) > self.maxTotFS.value:
				continue
			putEndSize = putEndSize + curUTSize
			grabInds.append(pi)
			if len(grabInds) >= maxOneGo:
				break
		# get keys (key:iv needs to be unique, key is unique, iv=0 should be fine)
		newFKeys = [tapeDB.getFolderKey(utSet[cfi]).hex().upper() for cfi in grabInds]
		comIVec = "00000000000000000000000000000000"
		# start tar and encrypt
		tarFPaths = [os.path.join(self.tempRoot.value,utSet[cfi] + ".tar.gz.aes") for cfi in grabInds]
		tarFOFs = [open(ctf,"wb") for ctf in tarFPaths]
		spTars = []
		spGzips = []
		spOssls = []
		for i in range(len(grabInds)):
			pi = grabInds[i]
			curSrcF = utSet[pi]
			curDstF = tarFPaths[i]
			spTar = subprocess.Popen(["tar","cp",curSrcF],stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,cwd=self.foldRoot.value)
			spGzip = subprocess.Popen(["gzip"],stdin=spTar.stdout,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
			spOssl = subprocess.Popen(["openssl","enc","-aes-256-ctr","-K",newFKeys[i],"-iv",comIVec],stdin=spGzip.stdout,stdout=tarFOFs[i],stderr=subprocess.DEVNULL)
			spTars.append(spTar)
			spGzips.append(spGzip)
			spOssls.append(spOssl)
		# wait for them to finish
		rcTars = [csp.wait() for csp in spTars]
		rcGzips = [csp.wait() for csp in spGzips]
		rcOssls = [csp.wait() for csp in spOssls]
		# if any failed, error out
		for i in range(len(grabInds)):
			curDstF = tarFPaths[i]
			if (rcTars[i] != 0) or (rcGzips[i] != 0) or (rcOssls[i] != 0):
				raise IOError("Problem writing " + curDstF)
		# add the actual final size
		actGrabInds = []
		realEndSize = startSize
		for i in range(len(grabInds)):
			curTSize = os.path.getsize(tarFPaths[i])
			if (realEndSize + curTSize) > self.maxTotFS.value:
				continue
			realEndSize = realEndSize + curTSize
			actGrabInds.append(grabInds[i])
		# move things around (not good, might make faster)
		i = len(actGrabInds)
		while i > 0:
			i = i - 1
			pi = actGrabInds[i]
			ntSet.append(utSet[pi])
			utSet.pop(pi)
			utSize.pop(pi)
		return realEndSize


class UnpackTapeProgram(whodunargs.StandardProgram):
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "UnpackTape"
		self.summary = "Restore data from tape."
		self.usage = "python3 TapeDeck UnpackTape"
		self.version = commonLicense
		self.dataRoot = whodunargs.ArgumentOptionFolderRead("--database","Database Directory","The folder containing the database.","--database C:\\Database\\")
		self.foldRoot = whodunargs.ArgumentOptionFolderWrite("--backfold","Backup Folder","The folder to back up.","--backfold C:\\WorkingSet\\")
		self.tempRoot = whodunargs.ArgumentOptionFolderWrite("--scratch","Scratch Directory","A folder for temporary tar files.","--scratch /tmp/vhs")
		self.reportFN = whodunargs.ArgumentOptionFileRead("--in","Input","The place to read tar backup.","--in /dev/st0")
		self.options.append(self.dataRoot)
		self.options.append(self.foldRoot)
		self.options.append(self.tempRoot)
		self.options.append(self.reportFN)
	def baseRun(self):
		if len(self.reportFN.value.strip()) == 0:
			raise ValueError("UnpackTape cannout read from stdin.")
		# load the database
		tapeDB = tapecom.TapeDatabase(self.dataRoot.value)
		knowFs = set(tapeDB.getFolders())
		# and the staging ground
		stageF = tapecom.StagingFolder(self.foldRoot.value)
		# make the scratch
		os.makedirs(self.tempRoot.value, exist_ok=False)
		# untar
		spTar = subprocess.Popen(["tar","xf",os.path.abspath(self.reportFN.value)],stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,cwd=self.tempRoot.value)
		if spTar.wait() != 0:
			raise IOError("Problem unpacking " + self.reportFN.value)
		# figure out which folders are present (complain if anything weird)
		upTars = os.listdir(self.tempRoot.value)
		upTarPaths = [os.path.join(self.tempRoot.value,cfn) for cfn in upTars]
		upFolds = []
		for i in range(len(upTars)):
			if not os.path.isfile(upTarPaths[i]):
				raise IOError(upTarPaths[i] + " is not a file, you sure you got the right tape?")
			if not upTars[i].endswith(".tar.gz.aes"):
				raise IOError(upTarPaths[i] + " extension does not match expected.")
			curFold = upTars[i][0:-11]
			if not (curFold in knowFs):
				raise IOError(upTarPaths[i] + " does not match a known folder.")
			upFolds.append(curFold)
		# of those present, figure out which are not already in staging
		stageFs = set(stageF.hotFolds)
		cpTarSrcs = []
		cpTarDsts = []
		for i in range(len(upTars)):
			if upFolds[i] in stageFs:
				continue
			cpTarSrcs.append(upTarPaths[i])
			cpTarDsts.append(upFolds[i])
		# unpack to staging
		maxOneGo = multiprocessing.cpu_count()
		while len(cpTarSrcs) > 0:
			# grab some stuff
			curUTSet = []
			curUTDst = []
			while len(curUTSet) < maxOneGo:
				if len(cpTarSrcs) == 0:
					break
				curUTSet.append(cpTarSrcs.pop(-1))
				curUTDst.append(cpTarDsts.pop(-1))
			# grab the keys
			newFKeys = [tapeDB.getFolderKey(cv).hex().upper() for cv in curUTDst]
			comIVec = "00000000000000000000000000000000"
			# set up jobs to untar
			spOssls = []
			spGzips = []
			spTars = []
			for i in range(len(curUTSet)):
				curSrcF = curUTSet[i]
				curDstF = curUTDst[i]
				spOssl = subprocess.Popen(["openssl","enc","-d","-aes-256-ctr","-K",newFKeys[i],"-iv",comIVec,"-in",curSrcF],stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
				spGzip = subprocess.Popen(["gzip","-d"],stdin=spOssl.stdout,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
				spTar = subprocess.Popen(["tar","xp"],stdin=spGzip.stdout,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,cwd=self.foldRoot.value)
				spTars.append(spTar)
				spGzips.append(spGzip)
				spOssls.append(spOssl)
			# wait for them to finish
			rcTars = [csp.wait() for csp in spTars]
			rcGzips = [csp.wait() for csp in spGzips]
			rcOssls = [csp.wait() for csp in spOssls]
			# complain if problem (trash the folder if problem)
			probFold = None
			for i in range(len(curUTSet)):
				if (rcTars[i] == 0) and (rcGzips[i] == 0) and (rcOssls[i] == 0):
					continue
				probFold = os.path.join(self.foldRoot.value,curUTDst[i])
				if os.path.exists(probFold):
					shutil.rmtree(probFold)
			if not (probFold is None):
				raise IOError("Problem unpacking " + probFold)
		# delete scratch
		shutil.rmtree(self.tempRoot.value)


class BackupDatabaseProgram(whodunargs.StandardProgram):
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "BackupDatabase"
		self.summary = "Backup the database."
		self.usage = "python3 TapeDeck BackupDatabase"
		self.version = commonLicense
		self.dataRoot = whodunargs.ArgumentOptionFolderRead("--database","Database Directory","The folder containing the database.","--database C:\\Database\\")
		self.reportFN = whodunargs.ArgumentOptionFileWrite("--out","Output","The place to write the final tar file.","--out /home/ben/backup.tar.gz.aes")
		self.usePass = whodunargs.ArgumentOptionString("--pass","Password","The password to encrypt with.","--pass 1337")
		self.options.append(self.dataRoot)
		self.options.append(self.reportFN)
		self.options.append(self.usePass)
		self.usePass.value = "1337"
	def baseRun(self):
		if len(self.reportFN.value.strip()) == 0:
			raise ValueError("BackupDatabase cannout output data on stdout.")
		# list all the contents of the directory
		allTarFs = os.listdir(self.dataRoot.value)
		# set up the commands
		tarFOF = open(self.reportFN.value, "wb")
		spTar = subprocess.Popen(["tar","cp"] + allTarFs,stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,cwd=self.dataRoot.value)
		spGzip = subprocess.Popen(["gzip"],stdin=spTar.stdout,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
		spOssl = subprocess.Popen(["openssl","enc","-aes-256-ctr","-pass","pass:" + self.usePass.value,"-pbkdf2"],stdin=spGzip.stdout,stdout=tarFOF,stderr=subprocess.DEVNULL)
		# wait for them to finish
		rcTar = spTar.wait()
		rcGzip = spGzip.wait()
		rcOssl = spOssl.wait()
		if (rcTar != 0) or (rcGzip != 0) or (rcOssl != 0):
			raise IOError("Problem backing up.")


class RestoreDatabaseProgram(whodunargs.StandardProgram):
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "RestoreDatabase"
		self.summary = "Restore the database."
		self.usage = "python3 TapeDeck RestoreDatabase"
		self.version = commonLicense
		self.dataRoot = whodunargs.ArgumentOptionFolderWrite("--database","Database Directory","The folder to extract the database to.","--database C:\\Database\\")
		self.reportFN = whodunargs.ArgumentOptionFileRead("--in","Input","The tar to backup from.","--in /home/ben/backup.tar.gz.aes")
		self.usePass = whodunargs.ArgumentOptionString("--pass","Password","The password to decrypt with.","--pass 1337")
		self.options.append(self.dataRoot)
		self.options.append(self.reportFN)
		self.options.append(self.usePass)
		self.usePass.value = "1337"
	def baseRun(self):
		if len(self.reportFN.value.strip()) == 0:
			raise ValueError("BackupDatabase cannout output data on stdout.")
		# figure out the directory to work from
		dbParent = os.path.abspath(self.dataRoot.value)
		# make damn sure the parent exists
		os.makedirs(dbParent, exist_ok=False)
		# set up the commands
		spOssl = subprocess.Popen(["openssl","enc","-d","-aes-256-ctr","-pass","pass:" + self.usePass.value,"-pbkdf2","-in",self.reportFN.value],stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
		spGzip = subprocess.Popen(["gzip","-d"],stdin=spOssl.stdout,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
		spTar = subprocess.Popen(["tar","xp"],stdin=spGzip.stdout,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,cwd=dbParent)
		# wait for them to finish
		rcTar = spTar.wait()
		rcGzip = spGzip.wait()
		rcOssl = spOssl.wait()
		if (rcTar != 0) or (rcGzip != 0) or (rcOssl != 0):
			raise IOError("Problem unpacking.")


class KillTapedProgram(whodunargs.StandardProgram):
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "KillTaped"
		self.summary = "Delete folders that are currently taped."
		self.usage = "python3 TapeDeck KillTaped"
		self.version = commonLicense
		self.dataRoot = whodunargs.ArgumentOptionFolderRead("--database","Database Directory","The folder containing the database.","--database C:\\Database\\")
		self.foldRoot = whodunargs.ArgumentOptionFolderWrite("--backfold","Backup Folder","The folder to clean up.","--backfold C:\\WorkingSet\\")
		self.reportFN = whodunargs.ArgumentOptionFileWrite("--out","Output","The place to list deleted folders.","--out /home/ben/report.txt")
		self.options.append(self.dataRoot)
		self.options.append(self.foldRoot)
		self.options.append(self.reportFN)
	def baseRun(self):
		# load the database
		tapeDB = tapecom.TapeDatabase(self.dataRoot.value)
		# get taped in staging folder
		stageF = tapecom.StagingFolder(self.foldRoot.value)
		utSet = stageF.getTaped(tapeDB)
		# open the output file
		outF = _openOutput(self.reportFN.value,self.useOut)
		# report
		for cf in utSet:
			shutil.rmtree(os.path.join(self.foldRoot.value,cf))
			outF.write(bytes(cf,"utf-8"))
			outF.write(bytes("\n","utf-8"))
		# close the file
		_closeOutput(self.reportFN.value,outF)


class FindTapeProgram(whodunargs.StandardProgram):
	def __init__(self):
		whodunargs.StandardProgram.__init__(self)
		self.name = "FindTape"
		self.summary = "Find the tape that holds a folder."
		self.usage = "python3 TapeDeck FindTapeListUntaped"
		self.version = commonLicense
		self.dataRoot = whodunargs.ArgumentOptionFolderRead("--database","Database Directory","The folder containing the database.","--database C:\\Database\\")
		self.reportFN = whodunargs.ArgumentOptionFileWrite("--out","Output","The place to list newly found folders.","--out /home/ben/report.tsv")
		self.limitFNs = whodunargs.ArgumentOptionStringGreedyVector("--fold","Limit Folders","The folders to look for. By default, all.","--fold a b")
		self.options.append(self.dataRoot)
		self.options.append(self.reportFN)
		self.options.append(self.limitFNs)
	def baseRun(self):
		# load the database
		tapeDB = tapecom.TapeDatabase(self.dataRoot.value)
		# get the tape inverse map
		invMap = tapeDB.getTapeInverseContentMap()
		# get the list of things to hunt
		huntFolds = self.limitFNs.value
		if len(huntFolds) == 0:
			huntFolds = tapeDB.getFolders()
		else:
			allFolds = set(tapeDB.getFolders())
			for cf in huntFolds:
				if not cf in allFolds:
					raise ValueError(cf + " is not a known folder.")
		# open output
		outF = _openOutput(self.reportFN.value,self.useOut)
		# output
		outF.write(bytes("Folder\tTape\n","utf-8"))
		for huntN in huntFolds:
			tapeN = "--UNASSIGNED--"
			if huntN in invMap:
				tapeN = invMap[huntN]
			outF.write(bytes(huntN,"utf-8"))
			outF.write(bytes("\t","utf-8"))
			outF.write(bytes(tapeN,"utf-8"))
			outF.write(bytes("\n","utf-8"))
		# close output
		_closeOutput(self.reportFN.value,outF)



class TapeDeckProgramSet(whodunargs.StandardProgramSet):
	def __init__(self):
		whodunargs.StandardProgramSet.__init__(self)
		self.name = "TapeDeck"
		self.summary = "Utilities for managing a backup."
		self.version = commonLicense
		self.programs["InitDatabase"] = lambda: InitDatabaseProgram()
		self.programs["ListFolders"] = lambda: ListFoldersProgram()
		self.programs["HuntNewFolders"] = lambda: HuntNewProgram()
		self.programs["VerifyFolders"] = lambda: VerifyFoldersProgram()
		self.programs["ListTapes"] = lambda: ListTapesProgram()
		self.programs["ListUntaped"] = lambda: ListUntapedProgram()
		self.programs["NewTape"] = lambda: NewTapeProgram()
		self.programs["UnpackTape"] = lambda: UnpackTapeProgram()
		self.programs["BackupDatabase"] = lambda: BackupDatabaseProgram()
		self.programs["RestoreDatabase"] = lambda: RestoreDatabaseProgram()
		self.programs["KillTaped"] = lambda: KillTapedProgram()
		self.programs["FindTape"] = lambda: FindTapeProgram()


if __name__ == "__main__":
	tps = TapeDeckProgramSet()
	wantR = tps.parseArguments(sys.argv[1:])
	if not (wantR is None):
		wantR.run()



'''
ListFolders
HashFolders
VerifyFolders
ListTapes
FindUntaped
NewTape
UnpackTape
BackupDatabase
RestoreDatabase
KillTaped
--
GetFolderTape

DropTape
'''


