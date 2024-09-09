PY := ${PINGVIN}/py-env/bin/python

prep-slice-grid:
	$(PY) generate-grid-2d.py 50 qcirc

prep-xsect :
	$(PY) generate-cross-sections.py 50 100

march :
	pvn march

vtk :
	pvn2vtk

clean:
	- rm qcirc.*
	- rm -rf xsect/
	- rm -rf pvnsim/


	

