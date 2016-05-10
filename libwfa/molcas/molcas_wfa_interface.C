#include <stdio.h>
#include <libwfa/libwfa.h>
extern "C"
{
    void molcas_wfa_interface_();
}

void molcas_wfa_interface_()
{
    std::cout << "C++ molcas_wfa_interface called" << std::endl;
}