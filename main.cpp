#include <iostream>
#include <gmp.h>
#include <algorithm>
#include <gmpxx.h>

void btrim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    [](unsigned char ch){return !std::isspace(ch);}));
}
void etrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         [](unsigned char ch) {return !std::isspace(ch);}).base(), s.end());
}
void trim(std::string& s) {
    etrim(s);
    btrim(s);
}
bool seq(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    return std::equal(a.begin(), a.end(), b.begin(), [](char a_char, char b_char) {
        return std::tolower(static_cast<unsigned char>(a_char)) ==
               std::tolower(static_cast<unsigned char>(b_char));
    });
}

int main(int argc, char* argv[]) {
    std::string input;
    std::cout << "/!\\ this might take a while..." << std::endl;
    bool mainloop = true;
    while (mainloop) {
        std::cout << "[E]ncode, [C]rack or [O]ther?: ";
        std::cin >> input;
        trim(input);
        if (seq(input, "e") || seq(input, "encode"))  {

            //set e
            mpz_class e("65537");
            std::cout << "Choose an exponent e (65537 if empty): ";
            std::cin >> input;
            trim(input);
            if (!seq(input, "")) e=input;
            else std::cout << "No input, defaulting to 65537..." << std::endl;
        }

        if (seq(input, "c") || seq(input, "crack")) {

        }

        if (seq(input, "o") || seq(input, "other")) {

        }
    }
}