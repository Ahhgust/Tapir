rm -rf test
python3 TapeDeck.py InitDatabase --database test/db
mkdir test/stage
mkdir test/stage/Thing1
echo hello > test/stage/Thing1/fileA.txt
python3 TapeDeck.py HuntNewFolders --database test/db --backfold test/stage
python3 TapeDeck.py ListUntaped --database test/db --backfold test/stage
python3 TapeDeck.py NewTape --database test/db --backfold test/stage --scratch test/vhs --out test/tape0.tar --name tape0
python3 TapeDeck.py KillTaped --database test/db --backfold test/stage/
python3 TapeDeck.py UnpackTape --database test/db --backfold test/stage --scratch test/vhs --in test/tape0.tar
python3 TapeDeck.py VerifyFolders --database test/db --backfold test/stage
python3 TapeDeck.py BackupDatabase --database test/db --out test/db_back.tar.gz.aes --pass 1337
rm -rf test/db
python3 TapeDeck.py RestoreDatabase --database test/db --in test/db_back.tar.gz.aes --pass 1337

