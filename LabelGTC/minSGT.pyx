from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "SuperGeneTrees/minSGT.h":
	string DoSuperGeneTree(string gcontent, string scontent, bool preserveDupSpec, string outputmode, string clades_to_preserve, string treated_trees)


cpdef getMinSGT(gcontent, scontent, bool preserveDupSpec, clades, trees, outputmode=""):
	
	cdef string outmode = outputmode



	print("----------------------------------------------------")
	res = DoSuperGeneTree(gcontent, scontent, preserveDupSpec, clades, trees, outmode)
	print("________________________________________________________")
	return res
