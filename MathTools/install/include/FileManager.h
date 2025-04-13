//This serves as a wrapper for c++ standard filesystem 
#pragma once
#include <filesystem>
#include <string>
void myCreateDirectory(std::string dir_name);

bool myFileExists(std::string file_name);