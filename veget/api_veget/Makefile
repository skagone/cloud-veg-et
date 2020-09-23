current:	delaware
hi:
	echo git reset HEAD~2   the last two commits
run:
	python3 api_veget.py 40N-100E
    
delaware:
	python3 api_veget.py -c delaware_config -s DRB.shp lilDRB_temp

test_cache:
	time python3 api_veget.py -c grid_config -s ./AOI/tile40N-110E_chip40N-102E.shp  tile40N-110E_chip40N-102E  2>&1 | tee  ./log/tile40N-110E_chip40N-102E.log&

help:
	python3 api_veget.py --help

git:
	echo Call Steffi and have her run make git in the root dir for this repo
	git remote set-url origin git@github.com:skagone/cloud-veg-et.git
	git config --global user.email skagone@contractor.usgs.gov
	git config --global user.name skagone
	git config --global push.default simple
	git add .
	git commit -m "automatic git update from Makefile"
	git push

graph:
	pyan3 --dot api_veget.py | dot -Tsvg -o api.svg
	pyan3 --dot ../vegetLib/vegetLib/veget.py | dot -Tsvg -o veget.svg

readme:
	echo https://github.com/tonybutzer/training/blob/master/craft/tools/README.md
