CC = gcc
BIN = ~/DGC_software/bin
LIBS	= -lm -lc 
C_OPTS	= -g

spillover:
%:	
	$(CC) $*.c -o $(BIN)/$@ $(C_OPTS) $(LIBS) 
	chmod  a+x  $(BIN)/$@


upload_version_starting_from_scratch:
	#(git init; git remote add spillover https://github.com/danigeos/spillover; git add .; git commit -a -mnewVersion; git push -u -f spillover master)


upload:
	#for initialization:  
	#git init; git remote add tisc https://github.com/danigeos/tisc; git add Makefile README config.mk bin demo doc include lib script src
	(git commit -a ; git push -u -f spillover master)


