#include <iostream>
#include <stdlib.h>
#include "H5Cpp.h"

using namespace H5;

int main(int argc, char** argv)
{
    std::cout << "Starting analysis of Molcas HDF5 file..." << std::endl;
    
    H5std_string file_name("molcas.rasscf.h5");
    if (argc>=2){
        file_name = argv[1];
    }
    std::cout << "Analysing file " << file_name << std::endl;
    
    H5File file( file_name, H5F_ACC_RDONLY ); // read-only for now
}