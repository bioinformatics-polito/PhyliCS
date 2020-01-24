CC=g++ 

cmd=$(CC) -o scripts/$@ scripts/$@.cpp

all: testBED status interval binUnsorted CNVcaller

testBED: scripts/testBED.cpp
	$(cmd)

status: scripts/status.cpp
	$(cmd)

interval: scripts/interval.cpp
	$(cmd)

binUnsorted: scripts/binUnsorted.cpp
	$(cmd)

CNVcaller: scripts/CNVcaller.cpp
	$(cmd)

