PY := ${PINGVIN}/py-env/bin/python

prep-slice-grid:
	$(PY) generate-grid-2d.py 50 qcirc

prep-xsect :
	$(PY) generate-cross-sections.py 50 100

march :
	pvn march

vtk :
	pvn2vtk

exit-plane:
	pvn2oned -@ 6.44

one-dimensionalise:
	pvn2oned --n-slices 50

plots:
	python plot-oned-props.py
	gnuplot plot-comparison.gplot

clean:
	- rm qcirc.*
	- rm -rf xsect/
	- rm -rf pvnsim/
