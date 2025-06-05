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
        std::getline(std::cin, input);
        trim(input);
        if (seq(input, "e") || seq(input, "encode"))  {

            //set e
            mpz_class e("65537");
            std::cout << "Choose an exponent e (65537 if empty): ";
            std::getline(std::cin, input);
            trim(input);
            if (!seq(input, "")) e=input;
            else std::cout << "No input, defaulting to 65537..." << std::endl;

            // setting p
            mpz_class p;
            while (true) {
                std::cout << "Choose p, a prime number: ";
                std::getline(std::cin, input);
                trim(input);
                p = input;
                // Wrong usage of method
                //TODO: fix
                if (mpz_probab_prime_p(p.get_mpz_t(), 30)==1||mpz_probab_prime_p(p.get_mpz_t(), 30)==2) {
                    p=input;
                    break;
                }
                else std::cout << "Not a prime number, try again..." << std::endl;
            }

            // setting q
            mpz_class q;
            while (true) {
                std::cout << "Choose q, a prime number: ";
                std::getline(std::cin, input);
                trim(input);
                q = input;
                // Wrong usage of method
                // TODO:fix
                if (mpz_probab_prime_p(q.get_mpz_t(), 30)==1||mpz_probab_prime_p(q.get_mpz_t(), 30)==2) {
                    q=input;
                    break;
                }
                else std::cout << "Not a prime number, try again..." << std::endl;
            }
            // set n
            mpz_class n=p*q;
            // set phi(n)
            mpz_class phi=(p-1)*(q-1);
            mpz_class d;
            mpz_invert(d.get_mpz_t(), e.get_mpz_t(), phi.get_mpz_t());
            std::cout << "Public key: (e = " << e << ", n = " << n << ")" << std::endl;
            std::cout << "Private key: (d = " << d << ", n = " << n << std::endl;
            while (true) {
                std::cout << "Enter a message to encrypt: ";
                std::getline(std::cin, input);
                trim(input);
                if (input.empty()) {
                    std::cout << "No input, exiting..." << std::endl;
                    break;
                }
                else {
                    // TODO: inplement input of strings
                    mpz_class m(input);
                    mpz_class c;
                    mpz_powm(c.get_mpz_t(), m.get_mpz_t(), e.get_mpz_t(), n.get_mpz_t());
                    std::cout << "Encrypted message: " << c << std::endl;
                    break;
                }
            }
        }

        if (seq(input, "c") || seq(input, "crack")) {

            //set e
            mpz_class e("65537");
            std::cout << "Choose an exponent e (65537 if empty): ";
            std::getline(std::cin, input);
            trim(input);
            if (!seq(input, "")) e=input;
            else std::cout << "No input, defaulting to 65537..." << std::endl;

            // set n
            mpz_class n;
            std::cout << "Enter n from the public key: ";
            trim(input);
            std::getline(std::cin, input);
            n = input;
            std::cout << "Trying to facorize n, this might take a while..." << std::endl;


            // TODO: multithread this
            //factorize n
            mpz_class p("2");
            mpz_class q("2");
            mpz_class max;
            mpz_class next;
            mpz_class nmod(n);
            mpz_sqrt(max.get_mpz_t(), n.get_mpz_t());
            while (p<=max){
                //try to find q
                if (mpz_divisible_p(n.get_mpz_t(), p.get_mpz_t()))  {
                    q = n/p;
                    if (mpz_probab_prime_p(q.get_mpz_t(), 30)==1||mpz_probab_prime_p(q.get_mpz_t(), 30)==2) {
                        std::cout << "Found p and q!" << std::endl;
                        std::cout << "p = " << p << std::endl;
                        std::cout << "q = " << q << std::endl;
                        break;
                    }
                }
                //get next prime
                mpz_nextprime(next.get_mpz_t(), p.get_mpz_t());
                p = next;
            }
            mpz_class phi((p-1)*(q-1));
            mpz_class d;
            mpz_invert(d.get_mpz_t(), e.get_mpz_t(), phi.get_mpz_t());
            std::cout << "Public key: (e = " << e << ", n = " << n << ")" << std::endl;
            std::cout << "Private key: (d = " << d << ", n = " << n << std::endl;
            while (true) {
                std::cout << "Enter a message to decrypt: ";
                std::getline(std::cin, input);
                trim(input);
                mpz_class c(input);
                mpz_class m;
                mpz_powm(m.get_mpz_t(), c.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());

                if (input.empty()) {
                    break;
                }
            }

        }

        if (seq(input, "o") || seq(input, "other")) {
            bool otherLoop= true;
            while (otherLoop) {
                std::cout << "Choose an operation: " << std::endl;
                std::cout << "     [1] Big number calculator" << std::endl;
                std::cout << "     [2] String to int (ASCII)" << std::endl;
                std::cout << "     [3] int (ASCII) to String" << std::endl;
                std::cout << "     [Q] Quit" << std::endl;
                std::cout << "Enter your choice: ";
                std::getline(std::cin, input);
                trim(input);
                switch (input[0]) {
                    case '1': {
                        std::cout << "Enter num1: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class num1(input);
                        std::cout << "Enter num2: ";
                        std::getline(std::cin, input);
                        mpz_class num2(input);
                        std::cout << "Choose an operation (+, -, *, /): ";
                        std::getline(std::cin, input);
                        trim(input);
                        switch (input[0]) {
                            case '+':
                                std::cout << num1 + num2 << std::endl;
                                break;
                            case '-':
                                std::cout << num1 - num2 << std::endl;
                                break;
                            case '*':
                                std::cout << num1 * num2 << std::endl;
                                break;
                            case '/':
                                std::cout << num1 / num2 << std::endl;
                                break;
                            default:
                                std::cout << "Invalid operation" << std::endl;
                                break;
                        }
                        otherLoop = false;
                        break;
                    }
                    case '2':
                        break;
                    case '3':
                        break;
                    case 'q':
                        mainloop = false;
                        break;
                    case 'Q':
                        mainloop = false;
                        break;
                    default:
                        std::cout << "Invalid choice" << std::endl;
                        break;
                }
            }
            continue;
        }
        else {
            std::cout << "Invalid choice" << std::endl;
            continue;
        }
    }
}