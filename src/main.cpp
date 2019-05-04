#include <iostream>

int globvar;

int main( int argc, char** argv ){
    globvar = 2;
    std::cout << globvar << std::endl;
    return 0;
}
