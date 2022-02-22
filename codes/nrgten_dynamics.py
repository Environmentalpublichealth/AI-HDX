"""
Created on August 2nd, 2021
Author: Jiali
Usage: compute dynamical signature of a protein based on protein structure
python nrgten_dynamics.py <input.pdb>
"""
import sys
input_file = sys.argv[1]

# need to install nrgten in the environmenet first
from nrgten.encom import ENCoM

if __name__ == "__main__":
	model = ENCoM(input_file)
	dyna_sig = model.compute_bfactors()
	print(dyna_sig)
