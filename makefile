install:
	- mkdir inst
	odin build src -out:inst/pvn -o:speed
	cp src/pvn2vtk.sh inst/pvn2vtk
	cp src/pvn2oned.sh inst/pvn2oned
	cp src/pvn-compute-norms.sh inst/pvn-compute-norms
	chmod +x inst/pvn2vtk inst/pvn2oned inst/pvn-compute-norms
	cp src/python-programs/pvn2vtk.py inst/
	cp src/python-programs/pvn2oned.py inst/
	cp src/python-programs/pvn-compute-norms.py inst/

clean:
	- rm -rf inst/


