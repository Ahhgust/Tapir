'''Common crap for tapes.'''

import os
import time
import struct
import hashlib


'''
A tape database is a set of files.
tapes.bin - Size, and add-time of each tape.
folders.bin - Size and add-time of each folder.
keys.bin - Keys for each folder, if present.
hashes.bin - Hash data for each folder.
map.bin - Which folders are in which tape.

*.ind

A staging folder is just a folder of folders.


'''

class IndexedFileDir:
	'''Contain the contents of an indexed file.'''
	def __init__(self,ipath,dpath):
		'''
		Set up a new indexed file.
		@param ipath: The path to the index file.
		@param dpath: The path to the data file.
		'''
		self.ipath = ipath;
		'''The path to the index file.'''
		self.dpath = dpath;
		'''The path to the data file.'''
		self.contOffs = {}
		'''The offsets (and lengths) of each content.'''
		# make an empty if it does not exist
		if not os.path.exists(ipath):
			makeF = open(ipath, "wb")
			makeF.close()
		if not os.path.exists(dpath):
			makeF = open(dpath, "wb")
			makeF.close()
		# load the contents
		loadI = open(ipath, "rb")
		while True:
			curNL = loadI.read(8)
			if len(curNL) != 8:
				if len(curNL) != 0:
					raise IOError(ipath + ": Index truncated.")
				break
			curNL = struct.unpack(">Q",curNL)[0]
			curN = loadI.read(curNL)
			if len(curN) != curNL:
				raise IOError(ipath + ": Index truncated.")
			curN = str(curN, "utf-8")
			curOLB = loadI.read(16)
			if len(curOLB) != 16:
				raise IOError(ipath + ": Index truncated.")
			self.contOffs[curN] = struct.unpack(">QQ",curOLB)
		loadI.close()
	def loadFile(self,fname):
		'''
		Load a file by name.
		@param fname: The name of the file to load.
		@return: The bytes of the file.
		'''
		curCD = self.contOffs[fname]
		dataI = open(self.dpath,"rb")
		dataI.seek(curCD[0])
		loadD = dataI.read(curCD[1])
		dataI.close()
		if len(loadD) != curCD[1]:
			raise IOError(self.dpath + ": File truncated")
		return loadD
	def addFile(self,fname,fbytes):
		'''
		Add a file to this thing.
		@param fname: The name of the file.
		@param fbytes: The bytes of the file.
		'''
		curOff = os.path.getsize(self.dpath);
		curN = bytes(fname,"utf-8")
		dataO = open(self.dpath,"ab")
		dataO.write(fbytes)
		dataO.close()
		indxO = open(self.ipath,"ab")
		indxO.write(struct.pack(">Q",len(curN)))
		indxO.write(curN)
		indxO.write(struct.pack(">QQ",curOff,len(fbytes)))
		indxO.close()
		self.contOffs[fname] = (curOff,len(fbytes))


class DirectoryTreeHash:
	'''A hash of the contents of a directory tree.'''
	def __init__(self):
		'''Set up an empty.'''
		self.folds = set()
		'''The folders in the thing.'''
		self.files = {}
		'''The files in the thing, with their hashes.'''
		self.weirds = set()
		'''All the weird files in the thing.'''
	def hashPath(self,rootPath,curPath):
		'''
		Add the hashes of the contents of a path.
		@param rootPath: The path to start from.
		@param curPath: The current sub-path.
		'''
		fullPath = os.path.join(rootPath,curPath)
		curConts = os.listdir(fullPath)
		for sub in curConts:
			curSCP = os.path.join(curPath,sub)
			subFP = os.path.join(fullPath,sub)
			if os.path.isfile(subFP):
				curH = hashlib.sha256()
				curF = open(subFP,"rb")
				curB = curF.read(0x010000)
				while len(curB) > 0:
					curH.update(curB)
					curB = curF.read(0x010000)
				curF.close()
				self.files[curSCP] = curH.digest()
			elif os.path.isdir(subFP):
				self.folds.add(curSCP)
				self.hashPath(rootPath,curSCP)
			else:
				self.weirds.add(curSCP)
	def packHash(self):
		'''
		Pack the contents of this hash to bytes ready to write.
		@return: The bytes to write.
		'''
		toR = bytearray()
		# output folders
		toR.extend(struct.pack(">Q",len(self.folds)))
		for cf in self.folds:
			cfb = bytes(cf,"utf-8")
			toR.extend(struct.pack(">Q",len(cfb)))
			toR.extend(cfb)
		# output files
		toR.extend(struct.pack(">Q",len(self.files)))
		for cf in self.files:
			cfh = self.files[cf]
			cfb = bytes(cf,"utf-8")
			toR.extend(struct.pack(">Q",len(cfb)))
			toR.extend(cfb)
			toR.extend(cfh)
		# sero marker
		toR.extend(struct.pack(">Q",0))
		return bytes(toR)
	def unpackHash(self,fromB):
		'''
		Unpack the contents of a hash.
		@param fromB: The previously packed contents.
		'''
		ci = 0
		# get folders
		numFold = struct.unpack(">Q",fromB[ci:ci+8])[0]
		ci += 8
		for i in range(numFold):
			foldNL = struct.unpack(">Q",fromB[ci:ci+8])[0]
			ci += 8
			foldN = fromB[ci:ci+foldNL]
			ci += foldNL
			foldN = str(foldN, "utf-8")
			self.folds.add(foldN)
		# hash length
		hashL = hashlib.sha256().digest_size
		# get files
		numFile = struct.unpack(">Q",fromB[ci:ci+8])[0]
		ci += 8
		for i in range(numFile):
			fileNL = struct.unpack(">Q",fromB[ci:ci+8])[0]
			ci += 8
			fileN = fromB[ci:ci+fileNL]
			ci += fileNL
			fileH = fromB[ci:ci+hashL]
			ci += hashL
			fileN = str(fileN, "utf-8")
			self.files[fileN] = fileH
		# make sure the finale is zero
		if struct.unpack(">Q",fromB[ci:ci+8])[0] != 0:
			raise IOError("Invalid marker at end of hash pack.")
		ci += 8
		# and make sure we are at the finale
		if ci != len(fromB):
			raise IOError("Trailing junk after hash pack.")
	def getDiffs(self,othHash):
		'''
		Get a list of differences between two hashes.
		@param othHash: The other hash.
		@return: A list of tuples (path,prob).
		'''
		allDiffs = []
		# folders first
		foldMiss = self.folds.difference(othHash.folds)
		for cf in foldMiss:
			allDiffs.append((cf,"Folder Missing."))
		foldExt = othHash.folds.difference(self.folds)
		for cf in foldExt:
			allDiffs.append((cf,"Unknown Folder."))
		# then files
		selfFiles = set(self.files)
		othFiles = set(othHash.files)
		fileMiss = selfFiles.difference(othFiles)
		for cf in fileMiss:
			allDiffs.append((cf,"File Missing."))
		fileExt = othFiles.difference(selfFiles)
		for cf in fileExt:
			allDiffs.append((cf,"Unknown File."))
		comFiles = selfFiles.intersection(othFiles)
		for cf in comFiles:
			if self.files[cf] != othHash.files[cf]:
				allDiffs.append((cf,"Hash Mismatch"))
		return allDiffs

_passKeySize = 32

class TapeDatabase:
	'''A collection of tapes and things in tapes.'''
	def __init__(self,dbPath):
		'''
		Set up a database.
		@param dbPath: The path the database is (or should be) in.
		'''
		# figure the paths to the things of interest
		self.tapeFile = os.path.join(dbPath, "tapes.bin")
		'''The path to basic info on the tapes.'''
		self.tapeFileI = os.path.join(dbPath, "tapes.ind")
		'''Index path.'''
		self.foldFile = os.path.join(dbPath, "folders.bin")
		'''The path to basic info on the folders.'''
		self.foldFileI = os.path.join(dbPath, "folders.ind")
		'''Index path.'''
		self.keyFile = os.path.join(dbPath, "keys.bin")
		'''The path to keys for each folder.'''
		self.keyFileI = os.path.join(dbPath, "keys.ind")
		'''Index path.'''
		self.hashFile = os.path.join(dbPath, "hashes.bin")
		'''The path to hashes of each folder.'''
		self.hashFileI = os.path.join(dbPath, "hashes.ind")
		'''Index path.'''
		self.mapFile = os.path.join(dbPath, "map.bin")
		'''The path to what each tape holds.'''
		self.mapFileI = os.path.join(dbPath, "map.ind")
		'''Index path.'''
		# actually open them up
		if not (os.path.exists(dbPath)):
			os.makedirs(dbPath)
		self.tapeD = IndexedFileDir(self.tapeFileI,self.tapeFile)
		'''The tapes in play.'''
		self.foldD = IndexedFileDir(self.foldFileI,self.foldFile)
		'''The folders in play.'''
		self.keyD = IndexedFileDir(self.keyFileI,self.keyFile)
		'''The keys in play.'''
		self.hashD = IndexedFileDir(self.hashFileI,self.hashFile)
		'''The hashes in play.'''
		self.mapD = IndexedFileDir(self.mapFileI,self.mapFile)
		'''The tape assignments in play.'''
	def getFolders(self):
		'''
		Get a list of known folders in this thing.
		@return: List of folder names.
		'''
		return list(self.foldD.contOffs)
	def getFolderData(self,foldN):
		'''
		Get the data for a folder.
		@param foldN: The name of the folder.
		@return: The data.
		'''
		foldDD = self.foldD.loadFile(foldN)
		return self._unpackDataFile(foldDD)
	def getTapes(self):
		'''
		Get a list of tapes in this thing.
		@return: List of tape names.
		'''
		return list(self.tapeD.contOffs)
	def getTapeData(self,tapeN):
		'''
		Get the data for a tape.
		@param tapeN: The name of the tape.
		@return: The data.
		'''
		tapeDD = self.tapeD.loadFile(tapeN)
		return self._unpackDataFile(tapeDD)
	def getFolderKey(self,foldN):
		'''
		Get the key for a folder.
		@param foldN: The name of the folder.
		@return: The bytes of the key.
		'''
		if foldN in self.keyD.contOffs:
			return self.keyD.loadFile(foldN)
		if not (foldN in self.foldD.contOffs):
			raise ValueError(foldN + " is not a known folder.")
		newKey = os.urandom(_passKeySize)
		self.keyD.addFile(foldN,newKey)
		return newKey
	def getFolderHash(self,foldN):
		'''
		Get the hash data for a folder.
		@param foldN: The name of the folder.
		@return: The parsed hash data.
		'''
		foldHD = self.hashD.loadFile(foldN)
		toR = DirectoryTreeHash()
		toR.unpackHash(foldHD)
		return toR
	def getTapeContents(self,tapeN):
		'''
		Get the folders that make up a tape.
		@param tapeN: The name of the tape.
		@return: The names of the folders that make it up.
		'''
		allFold = []
		tapeCD = self.mapD.loadFile(tapeN)
		ci = 0
		while ci < len(tapeCD):
			curL = struct.unpack(">Q",tapeCD[ci:ci+8])[0]
			ci += 8
			curN = tapeCD[ci:ci+curL]
			ci += curL
			curN = str(curN, "utf-8")
			allFold.append(curN)
		return allFold
	def addFolder(self,foldN,foldDat,foldHash):
		'''
		Add a folder for consideration.
		@param foldN: The name of the folder.
		@param foldDat: The data map of the folder.
		@param foldHash: The hash data for the folder.
		'''
		if foldN in self.foldD.contOffs:
			raise IOError(foldN + " is already in database.")
		hBytes = foldHash.packHash()
		dBytes = self._packDataFile(foldDat)
		self.hashD.addFile(foldN,hBytes)
		self.foldD.addFile(foldN,dBytes)
	def addTape(self,tapeN,tapeDat,contFolds):
		'''
		Add a tape.
		@param tapeN: The name of the tape.
		@param tapeDat: The data map of the tape.
		@param contFolds: The folders in the tape.
		'''
		if tapeN in self.tapeD.contOffs:
			raise IOError(tapeN + " is already in database.")
		cBytes = bytearray()
		for cf in contFolds:
			cfb = bytes(cf,"utf-8")
			cBytes.extend(struct.pack(">Q",len(cfb)))
			cBytes.extend(cfb);
		dBytes = self._packDataFile(tapeDat)
		self.mapD.addFile(tapeN,bytes(cBytes))
		self.tapeD.addFile(tapeN,dBytes)
	def getTapeContentMap(self):
		'''
		Get the full map from tape to contents.
		@return: Map from tape name to contents.
		'''
		toR = {}
		for tapeN in self.tapeD.contOffs:
			toR[tapeN] = self.getTapeContents(tapeN)
		return toR
	def getTapeInverseContentMap(self):
		'''
		Get the full map from contents to tape.
		@return: Map from folder name to tape name.
		'''
		forwMap = self.getTapeContentMap()
		toR = {}
		for tapeN in forwMap:
			curFolds = forwMap[tapeN]
			for cf in curFolds:
				toR[cf] = tapeN
		return toR
	def _unpackDataFile(self,dataFC):
		'''
		Unpack the contents of a data file.
		@param dataFC: The bytes of the data file.
		@return: The unpacked map.
		'''
		toR = {}
		ci = 0
		while ci < len(dataFC):
			keyNL = struct.unpack(">Q",dataFC[ci:ci+8])[0]
			ci += 8
			keyN = dataFC[ci:ci+keyNL]
			ci += keyNL
			valNL = struct.unpack(">Q",dataFC[ci:ci+8])[0]
			ci += 8
			valN = dataFC[ci:ci+valNL]
			ci += valNL
			toR[str(keyN,"utf-8")] = str(valN,"utf-8")
		return toR
	def _packDataFile(self,dataFM):
		'''
		Pack up a map as a data file.
		@param dataFM: The map to pack up.
		@return: The packed bytes.
		'''
		toR = bytearray()
		for curN in dataFM:
			curK = bytes(curN,"utf-8")
			curV = bytes(dataFM[curN],"utf-8")
			toR.extend(struct.pack(">Q",len(curK)))
			toR.extend(curK)
			toR.extend(struct.pack(">Q",len(curV)))
			toR.extend(curV)
		return toR


def getFolderSize(foldP):
	'''
	Get the (approximate) size of a folder.
	@param foldP: The path to the folder.
	@return: The size in bytes.
	'''
	totalSize = 0
	for cf in os.listdir(foldP):
		curFP = os.path.join(foldP,cf)
		if os.path.isdir(curFP):
			totalSize = totalSize + getFolderSize(curFP)
		elif os.path.isfile(curFP):
			totalSize = totalSize + os.path.getsize(curFP)
	return totalSize


class StagingFolder:
	'''A folder things are put into before writing to tape.'''
	def __init__(self,rootPath):
		'''
		Set up a staging folder.
		@param rootPath: The path to the staging folder.
		'''
		self.rootPath = rootPath
		'''The path to the folder.'''
		self.hotFolds = os.listdir(rootPath)
		'''The folders in the path of interest.'''
		for cv in self.hotFolds:
			curP = os.path.join(rootPath, cv)
			if not os.path.isdir(curP):
				raise IOError(curP + " is not a folder.")
	def getUnknown(self,toTD):
		'''
		Note which folders in this thing are unknown to a tape database.
		@param toTD: The tape database to consider.
		@return: The set of unknown folders.
		'''
		knownFs = set(toTD.getFolders())
		hotFs = set(self.hotFolds)
		return hotFs.difference(knownFs)
	def getKnown(self,toTD):
		'''
		Note which folders in this thing are known to a tape database.
		@param toTD: The tape database to consider.
		@return: The list of known folders.
		'''
		knownFs = set(toTD.getFolders())
		hotFs = set(self.hotFolds)
		return hotFs.intersection(knownFs)
	def addUnknown(self,toTD):
		'''
		Add any folders in this to the tape database (if not already present).
		@param toTD: The tape database to add to.
		@return: The set of added folders.
		'''
		knownFs = set(toTD.getFolders())
		addedFs = set()
		for cf in self.hotFolds:
			if cf in knownFs:
				continue
			# get the hash of the thing
			cHash = DirectoryTreeHash()
			cHash.hashPath(self.rootPath,cf)
			# get the size of the thing
			totalSize = getFolderSize(os.path.join(self.rootPath,cf))
			# set up the map
			datMap = {}
			datMap["Size"] = repr(totalSize)
			datMap["Added"] = repr(int(time.time()))
			# add the thing
			toTD.addFolder(cf,datMap,cHash)
			addedFs.add(cf)
		return addedFs
	def getTaped(self,toTD):
		'''
		Get the folders in this thing that are currently in a tape.
		@param toTD: The tape database to consider.
		@return: The list of untaped folders.
		'''
		# get the set of things in tapes
		allTM = toTD.getTapeContentMap()
		allTS = []
		for tapeN in allTM:
			allTS.extend(allTM[tapeN])
		allTS = set(allTS)
		# get the things in this that are accounted for
		hotFs = set(self.hotFolds)
		return list(hotFs.intersection(allTS))
	def getUntaped(self,toTD):
		'''
		Get the folders in this thing that are not currently in a tape.
		@param toTD: The tape database to consider.
		@return: The list of untaped folders, in the order they were added to the database.
		'''
		# get the set of things in tapes
		allTM = toTD.getTapeContentMap()
		allTS = []
		for tapeN in allTM:
			allTS.extend(allTM[tapeN])
		allTS = set(allTS)
		# get the things in this that are unaccounted for
		hotFs = set(self.hotFolds)
		utFolds = list(hotFs.difference(allTS))
		# get the time each was added
		addFolds = []
		for uf in utFolds:
			if not (uf in toTD.foldD.contOffs):
				continue
			ufDat = toTD.getFolderData(uf)
			addT = int(ufDat["Added"])
			addFolds.append((addT,uf))
		# sort
		addFolds.sort()
		# and unpack
		return [cv[1] for cv in addFolds]
	def checkFolderHash(self,toTD,foldN):
		'''
		Check that the hash of a folder matches that in the database.
		@param toTD: The tape database to compare to.
		@param foldN: The name of the folder to look at.
		@return: The difference list.
		'''
		# get the known hash
		dataH = toTD.getFolderHash(foldN)
		# get the hash of the current crap
		fileH = DirectoryTreeHash()
		fileH.hashPath(self.rootPath,foldN)
		# return the differences
		return dataH.getDiffs(fileH)





