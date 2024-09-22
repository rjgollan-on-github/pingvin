install:
	- mkdir inst
	odin build src -out:inst/pvn -o:speed
	cp src/pvn2vtk.sh inst/pvn2vtk
	cp src/pvn2oned.sh inst/pvn2oned
	chmod +x inst/pvn2vtk inst/pvn2oned
	cp src/python-programs/pvn2vtk.py inst/
	cp src/python-programs/pvn2oned.py inst/

clean:
	- rm -rf inst/


