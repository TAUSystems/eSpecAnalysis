#pragma once

#ifndef __main_h__
#define __main_h__

#include <string>

// extern "C" is necessary so that the analyze_espec_images function name 
// does not get mangled in the shared library
#ifdef __cplusplus
extern "C" {
#endif

void analyze_espec_images(std::string filename); 

#ifdef __cplusplus
}
#endif

#endif
