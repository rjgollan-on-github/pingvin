install:
	- mkdir inst
	odin build src -out:inst/pvn -o:speed
	cp src/pvn2vtk.sh inst/pvn2vtk
	chmod +x inst/pvn2vtk
	cp src/python-programs/pvn2vtk.py inst/

clean:
	- rm -rf inst/


