#include"SihatApplication.h"
#include<iostream>
#include<stdexcept>

// --- Main Application ---
int main()
{
    try
    {
        SihatApplication app("../source/configFiles/default.json");
        app.run();
    }
    catch (const std::exception& e){
        std::cerr << "A fatal error occured: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}