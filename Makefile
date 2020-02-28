publish:
	git config --global user.email skagone@contractor.usgs.gov
	git config --global user.name skagone
	git config --global push.default simple
	git add .
	git commit -m "automatic git update from Makefile"
	git push


up:
	(cd swarm; make start; make xrdp)


lite:
	(cd opt; git clone http://github.com/tonybutzer/lite-stac)
