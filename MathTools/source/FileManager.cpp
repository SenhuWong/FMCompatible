#include "FileManager.h"
#include <iostream>
void myCreateDirectory(std::string dir_name)
{
    std::filesystem::path current_path = std::filesystem::current_path();
    std::filesystem::path dir_path = current_path/dir_name;
    bool alreadyExist = std::filesystem::is_directory(std::filesystem::status(dir_path));
    if(alreadyExist)
    {

    }
    else
    {
        std::filesystem::create_directory(dir_path);
    }
}

bool myFileExists(std::string file_name)
{
    std::filesystem::path current_path = std::filesystem::current_path();
    std::filesystem::path file_path = current_path/file_name;
    bool alreadyExist = std::filesystem::is_regular_file(std::filesystem::status(file_path));
    return alreadyExist;
}