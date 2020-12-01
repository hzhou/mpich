##
## Copyright (C) by Argonne National Laboratory
##     See COPYRIGHT in top-level directory
##

from local_python import MPI_API_Global as G
from local_python.mpi_api import load_mpi_api
from local_python.binding_c import dump_mpi_c, dump_Makefile_mk
from local_python import RE
import os
import glob

def main():
    print("Loading maint/mpi_standard_api.txt ...")
    load_mpi_api("maint/mpi_standard_api.txt")

    os.chdir("src/binding/c")
    api_files = glob.glob("*_api.txt")
    for f in api_files:
        if RE.match(r'(\w+)_api.txt', f):
            # The name in eg pt2pt_api.txt indicates the output folder.
            # Only the api functions with output folder will get generated.
            # This allows simple control of what functions to generate.
            print("Loading src/binding/c/%s ..." % f)
            load_mpi_api(f, RE.m.group(1))

    func_list = [f for f in G.FUNCS.values() if 'dir' in f]
    func_list.sort(key = lambda f: f['dir'])
    for func in func_list:
        if RE.search(r'not_implemented', func['attrs']):
            print("  skip %s (not_implemented)" % func['name'])
            pass
        else:
            dump_mpi_c(func)

    dump_Makefile_mk()

# ---------------------------------------------------------
if __name__ == "__main__":
    main()
