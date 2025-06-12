#include <iostream>
#include <gmp.h>
#include <algorithm>
#include <gmpxx.h>
#include <ctime>
#include <thread>
#include <vector>
#include <cmath>
#include <iomanip>
#include <bits/random.h>

#include "MontgomeryCurve.h"

//Choose-able prime: 100200300400500600700800900000000009008007006005004003002051

// Globals for thread communication
std::atomic<bool> found(false);
std::mutex result_mutex;
mpz_class final_p;
mpz_class final_q;
std::atomic<uint64_t> primes_checked(0);
    // Elliptic Curves
std::mutex cout_mutex;
std::atomic<bool> found_factor(false);
std::atomic<unsigned long long> total_curves(0);
std::mutex factor_mutex;

void factor_thread(mpz_class n, const mpz_class& start, const mpz_class& end) {
    mpz_class p;

    mpz_class start_minus1 = start - 1;
    mpz_nextprime(p.get_mpz_t(), start_minus1.get_mpz_t());

    while (p <= end && !found.load()) {
        ++primes_checked;
        if (mpz_divisible_p(n.get_mpz_t(), p.get_mpz_t())) {
            mpz_class q = n / p;
            if (mpz_probab_prime_p(q.get_mpz_t(), 30) >= 1) {
                std::lock_guard<std::mutex> lock(result_mutex);
                if(!found) {
                    final_p = p;
                    final_q = q;
                    found = true;
                }
                return;
            }
        }
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
    }
}

void ecm_thread(const mpz_class &n, const mpz_class &k_B1, const mpz_class &k_B2,
                const std::vector<mpz_class> &primes, gmp_randstate_t state, unsigned thread_id) {

    while (!found_factor.load()) {
        total_curves++;

        mpz_class A, x0;
        mpz_urandomm(x0.get_mpz_t(), state, n.get_mpz_t());
        mpz_urandomm(A.get_mpz_t(), state, n.get_mpz_t());
        mpz_class discriminant = (A*A - 4);
        mpz_mod(discriminant.get_mpz_t(), discriminant.get_mpz_t(), n.get_mpz_t());
        if (discriminant == 0) continue;

        MontgomeryCurve curve(A, n);
        MontgomeryPoint P(x0, 1);

        MontgomeryPoint result = curve.scalar_multiply(k_B1, P);
        mpz_class gcd;
        mpz_gcd(gcd.get_mpz_t(), result.Z.get_mpz_t(), n.get_mpz_t());

        if (gcd != 1 && gcd != n) {
            found_factor = true;
            {
                std::lock_guard<std::mutex> lock(factor_mutex);
                final_p = gcd;
                final_q = n / gcd;
            }
            {
                std::lock_guard<std::mutex> lock(cout_mutex);
                std::cout << "\nThread " << thread_id << ": Factors found in Phase 1 after " << total_curves.load() << " curves:\n";
                std::cout << "Factor p: " << final_p << "\n";
                std::cout << "Factor q: " << final_q << "\n";
            }
            return;
        }

        // Phase 2
        MontgomeryPoint Q = result;
        mpz_class gcd2;

        for (const auto &p : primes) {
            if (p <= k_B1) continue;
            if (p > k_B2) break;
            Q = curve.scalar_multiply(p, Q);
            mpz_gcd(gcd2.get_mpz_t(), Q.Z.get_mpz_t(), n.get_mpz_t());
            if (gcd2 != 1 && gcd2 != n) break;
        }

        if (gcd2 != 1 && gcd2 != n) {
            found_factor = true;
            {
                std::lock_guard<std::mutex> lock(factor_mutex);
                final_p = gcd2;
                final_q = n / gcd2;
            }
            {
                std::lock_guard<std::mutex> lock(cout_mutex);
                std::cout << "\nThread " << thread_id << ": Factors found in Phase 2 after " << total_curves.load() << " curves:\n";
                std::cout << "Factor p: " << final_p << "\n";
                std::cout << "Factor q: " << final_q << "\n";
            }
            return;
        }

        if (total_curves % 1000 == 0 && thread_id == 0) {
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::cout << "Thread " << thread_id << ": " << total_curves.load() << " curves tested...\n";
        }
    }
}


void progress_display(const size_t total_primes) {
    constexpr int barWidth = 50; // Width of the progress bar
    while (!found.load()) {
        uint64_t checked = primes_checked.load();
        double progress = static_cast<double>(checked) / static_cast<double>(total_primes);
        int pos = static_cast<int>(barWidth * progress);

        std::cout << "\r|";
        for (int i = 0; i < barWidth; ++i) {
            if (i <= pos) std::cout << "█";
            else std::cout << " ";
        }

        std::cout << "| " << std::fixed << std::setprecision(2)
                  << (progress * 100.0) << "%    " << std::flush;

        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }

    // Print complete bar when done
    std::cout << "\r|";
    for (int i = 0; i < barWidth; ++i) std::cout << "█";
    std::cout << "| 100.00%" << std::endl;
}

void trimStart(std::string& s) {
    s.erase(s.begin(), std::ranges::find_if(s,
                                            [](unsigned char ch){return !std::isspace(ch);}));
}
void trimEnd(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         [](unsigned char ch) {return !std::isspace(ch);}).base(), s.end());
}
void trim(std::string& s) {
    trimStart(s);
    trimEnd(s);
}
bool seq(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    return std::equal(a.begin(), a.end(), b.begin(), [](char a_char, char b_char) {
        return std::tolower(static_cast<unsigned char>(a_char)) ==
               std::tolower(static_cast<unsigned char>(b_char));
    });
}
mpz_class next_palindrome(const mpz_class& n) {
    std::string s = n.get_str();
    const int len = static_cast<int>(s.length());
    std::string left = s.substr(0, (len + 1) / 2);
    std::string mirrored = left;
    if (len % 2 == 0)
        mirrored += std::string(left.rbegin(), left.rend());
    else
        mirrored += std::string(left.rbegin() + 1, left.rend());

    mpz_class pal(mirrored);

    if (pal > n)
        return pal;
    mpz_class left_num(left);
    left_num += 1;
    left = left_num.get_str();

    std::string new_pal = left;
    if (len % 2 == 0)
        new_pal += std::string(left.rbegin(), left.rend());
    else
        new_pal += std::string(left.rbegin() + 1, left.rend());

    return mpz_class(new_pal);
}
std::string mpz_to_ascii_string(const mpz_class& num) {
    std::string result;
    mpz_class temp = num;
    while (temp > 0) {
        char c = static_cast<char>(mpz_class(temp % 256).get_ui());
        result += c;
        temp /= 256;
    }
    std::ranges::reverse(result); // fix order
    return result;
}
mpz_class ascii_string_to_mpz(const std::string& str) {
    mpz_class result = 0;
    for (char c : str) {
        result *= 256;
        result += static_cast<unsigned char>(c);
    }
    return result;
}
std::string convert_base(const std::string& number, int from_base, int to_base) {
    if (from_base < 2 || from_base > 62 || to_base < 2 || to_base > 62) {
        throw std::invalid_argument("Base must be between 2 and 62.");
    }
    mpz_class value;
    // Set value from the string using the from_base
    if (mpz_set_str(value.get_mpz_t(), number.c_str(), from_base) != 0) {
        throw std::invalid_argument("Invalid number for base " + std::to_string(from_base));
    }

    // Convert the value to a string in the target base
    return value.get_str(to_base);
}
size_t estimate_total_primes(const mpz_class& max) {
    double max_d = mpz_get_d(max.get_mpz_t());
    return static_cast<size_t>(max_d / std::log(max_d));
}

int main() {
    std::string input;
    std::cout << "/!\\ this might take a while..." << std::endl;
    bool mainloop = true;
    while (mainloop) {
        std::cout << "[E]ncode, [C]rack, [L]Elliptic Curve Cracking or [O]ther?: ";
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
            bool encryptLoop = true;
            while (encryptLoop) {
                std::cout << "do you want to encrypt an int or a string?" << std::endl;
                std::cout << "     [1] int" << std::endl;
                std::cout << "     [2] string" << std::endl;
                std::cout << "Enter your choice: ";
                std::getline(std::cin, input);
                trim(input);
                switch (input[0]) {
                    case '1': {
                        std::cout << "Enter the number: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class m(input);
                        mpz_class c;
                        mpz_powm(c.get_mpz_t(), m.get_mpz_t(), e.get_mpz_t(), n.get_mpz_t());
                        std::cout << "Encrypted message: " << c << std::endl;
                        encryptLoop = false;
                        break;
                    }
                    case '2': {
                        std::cout << "Enter the string: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class m(ascii_string_to_mpz(input));
                        mpz_class c;
                        mpz_powm(c.get_mpz_t(), m.get_mpz_t(), e.get_mpz_t(), n.get_mpz_t());
                        std::cout << "Encrypted message: " << c << std::endl;
                        encryptLoop = false;
                        break;
                    }
                    default: {
                        std::cout << "Invalid choice" << std::endl;
                    }
                }
            }
            continue;
        }


        if (seq(input, "l") || seq(input, "elliptic")) {
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

            // get B1 based on size of n
            //TODO: rework
            unsigned long int B1;
            mpz_class digitsOfFactor;
            std::cout << "How many digits does the smallest factor have? Enter to skip and choose approximate value: ";
            std::getline(std::cin, input);
            trim(input);
            if (!seq(input, "")) {
                digitsOfFactor = input;
                std::cout << "set the length of a factor of n to be approximately " << digitsOfFactor << " digits" << std::endl;
            }
            else {
                // guess factor size to be ~sqrt n
                mpz_class sqrt_n;
                mpz_sqrt(sqrt_n.get_mpz_t(), n.get_mpz_t());
                unsigned long int digits = mpz_sizeinbase(sqrt_n.get_mpz_t(), 10);
                digitsOfFactor = digits;
                std::cout << "guessed the length of a factor n to be approximately " << digitsOfFactor << " digits" << std::endl;
            }
            // Basierend auf geschätzter Faktorbitlänge B1 auswählen
            if (digitsOfFactor >= 40) B1 = 50000000;
            else if (digitsOfFactor >= 35) B1 = 10000000;
            else if (digitsOfFactor >= 30) B1 = 3000000;
            else if (digitsOfFactor >= 25) B1 = 1000000;
            else if (digitsOfFactor >= 20) B1 = 250000;
            else if (digitsOfFactor >= 15) B1 = 25000;
            else if (digitsOfFactor >= 10) B1 = 2000;
            else return 500;
            std::cout << "Based on length of factor of n chose B1 to be: " << B1 << std::endl;

            unsigned long int B2(50 * B1);

            // Calculate all primes up to B2
            std::vector<bool> is_prime(B2+1, true);
            is_prime[0] = is_prime[1] = false;
            for (unsigned long int i = 2; i * i <= B2; ++i) {
                if (is_prime[i]) {
                    for (unsigned long int j = i * i; j <= B2; j += i) {
                        is_prime[j] = false;
                    }
                }
            }
            std::vector<mpz_class> primes;
            for (unsigned long int i = 2; i <= B2; ++i) {
                if (is_prime[i]) primes.emplace_back(i);
            }
            std::cout << "found "<< primes.size() << " primes in range 2 to B2" << std::endl;



            // calculate k_Bx (the scalar for the point multiplication)
            mpz_class k_B1(1);
            mpz_class k_B2(1);
            // k = \prod_p^B (p)^(round-down to next int(log_p(B))) wobei p stets prim
            for (mpz_class p : primes) {
                if (p>B1) break;
                double exponent = std::floor(std::log(B1)/std::log(p.get_d()));
                mpz_class max_pow;
                mpz_pow_ui(max_pow.get_mpz_t(), p.get_mpz_t(), static_cast<unsigned long long>(exponent));
                mpz_lcm(k_B1.get_mpz_t(), k_B1.get_mpz_t(), max_pow.get_mpz_t());
            }
            for (mpz_class p : primes) {
                double exponent = std::floor(std::log(B2)/std::log(p.get_d()));
                mpz_class max_pow;
                mpz_pow_ui(max_pow.get_mpz_t(), p.get_mpz_t(), static_cast<unsigned long long>(exponent));
                mpz_lcm(k_B2.get_mpz_t(), k_B2.get_mpz_t(), max_pow.get_mpz_t());
            }
            std::cout << "k_B1: " << k_B1 << std::endl;
            auto curveBeginning = std::chrono::high_resolution_clock::now();

            // get amount of threads
            unsigned int num_threads = std::thread::hardware_concurrency();
            if (num_threads == 0) {
                std::cout << "No threads detected, deafulting to 4" << std::endl;
                num_threads = 4;
            } else std::cout << "\nDetected " << num_threads << " Threads" << std::endl;
            std::vector<std::thread> threads;

            // Launch threads
            for (unsigned i = 0; i < num_threads; ++i) {
                threads.emplace_back([&, i]() {
                    gmp_randstate_t local_state;
                    gmp_randinit_mt(local_state);
                    // Seed using random_device or time + thread id
                    unsigned long seed = std::random_device{}() + i * 7919;
                    gmp_randseed_ui(local_state, seed);

                    ecm_thread(n, k_B1, k_B2, primes, local_state, i);

                    gmp_randclear(local_state);
                });
            }

            for (auto &t : threads) {
                if (t.joinable()) t.join();
            }

            auto curveEnding = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> CurveElapsed = curveEnding - curveBeginning;
            long elapsedS = static_cast<long>(CurveElapsed.count());
            long elapsedHours = elapsedS / 3600;
            long elapsedMinutes = (elapsedS % 3600) / 60;
            long elapsedSeconds = elapsedS % 60;
            if (elapsedHours > 0) std::cout << "Factorizing n took: " << elapsedHours << " Hours, " << elapsedMinutes << " Minutes and " << elapsedSeconds << " Seconds" << std::endl;
            else if (elapsedMinutes > 0) std::cout << "Factorizing n took: " << elapsedMinutes << " Minutes and " << elapsedSeconds << " Seconds" << std::endl;
            else std::cout << "Factorizing n took: " << elapsedSeconds << " Seconds" << std::endl;

            mpz_class phi((final_p-1)*(final_q-1));
            mpz_class d;
            mpz_invert(d.get_mpz_t(), e.get_mpz_t(), phi.get_mpz_t());
            std::cout << "Public key: (e = " << e << ", n = " << n << ")" << std::endl;
            std::cout << "Private key: (d = " << d << ", n = " << n << ")" << std::endl;
            bool crackLoop = true;
            while (crackLoop) {
                std::cout << "Do you want to get the decoded message as an int or a string?" << std::endl;
                std::cout << "     [1] int" << std::endl;
                std::cout << "     [2] string" << std::endl;
                std::cout << "Enter your choice: ";
                std::getline(std::cin, input);
                trim(input);
                switch (input[0]) {
                    case '1': {
                        std::cout << "Enter the encrypted message: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class c(input);
                        mpz_class m;
                        mpz_powm(m.get_mpz_t(), c.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());
                        std::cout << "Decrypted message: " << m << std::endl;
                        crackLoop = false;
                        break;
                    }
                    case '2': {
                        std::cout << "Enter the encrypted message: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class c(input);
                        mpz_class m;
                        mpz_powm(m.get_mpz_t(), c.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());
                        std::cout << "Decrypted message: " << mpz_to_ascii_string(m) << std::endl;
                        crackLoop = false;
                        break;
                    }
                    default: {
                        std::cout << "Invalid choice" << std::endl;
                    }
                }
            }
            continue;
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
            auto beginning = std::chrono::high_resolution_clock::now();
            std::cout << "Trying to factorize n, this might take a while..." << std::endl;
            mpz_class max;
            mpz_sqrt(max.get_mpz_t(), n.get_mpz_t());

            unsigned int NUM_THREADS = std::thread::hardware_concurrency();
            if (NUM_THREADS == 0) {
                NUM_THREADS = 4;
                std::cout << "couldn't detect amount of threads, using 4 instead" << std::endl;
            }
            else
                std::cout << "detected " << NUM_THREADS << " threads" << std::endl;
            std::vector<std::thread> threads;
            mpz_class current("2");
            mpz_class range_step;

            mpz_class max_minus_2 = max - 2;
            mpz_class chunk_size = max_minus_2 / NUM_THREADS;

            std::thread progress_thread(progress_display, estimate_total_primes(max));
            std::cout << std::endl;
            found = false;

            // Launch Threads
            for (int i = 0; i < NUM_THREADS; i++) {
                mpz_class start = 2 + i * chunk_size;
                mpz_class end = (i==NUM_THREADS-1) ? max : start + chunk_size;

                threads.emplace_back(factor_thread, n, start, end);
            }

            for (auto& t: threads) t.join();

            if (!found) {
                std::cout << "failed to factorize n" << std::endl;
                return 1;
            }

            std::this_thread::sleep_for(std::chrono::milliseconds(250));
            std::cout << "\nFound p and q!" << std::endl;
            std::cout << "p = " << final_p << std::endl;
            std::cout << "q = " << final_q << std::endl;

            progress_thread.detach();
            auto ending = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = ending - beginning;
            long elapsedS = static_cast<long>(elapsed.count());
            long elapsedHours = elapsedS / 3600;
            long elapsedMinutes = (elapsedS % 3600) / 60;
            long elapsedSeconds = elapsedS % 60;
            if (elapsedHours > 0) std::cout << "Factorizing n took: " << elapsedHours << " Hours, " << elapsedMinutes << " Minutes and " << elapsedSeconds << " Seconds" << std::endl;
            else if (elapsedMinutes > 0) std::cout << "Factorizing n took: " << elapsedMinutes << " Minutes and " << elapsedSeconds << " Seconds" << std::endl;
            else std::cout << "Factorizing n took: " << elapsedSeconds << " Seconds" << std::endl;

            mpz_class phi((final_p-1)*(final_q-1));
            mpz_class d;
            mpz_invert(d.get_mpz_t(), e.get_mpz_t(), phi.get_mpz_t());
            std::cout << "Public key: (e = " << e << ", n = " << n << ")" << std::endl;
            std::cout << "Private key: (d = " << d << ", n = " << n << std::endl;
            bool crackLoop = true;
            while (crackLoop) {
                std::cout << "Do you want to get the decoded message as an int or a string?" << std::endl;
                std::cout << "     [1] int" << std::endl;
                std::cout << "     [2] string" << std::endl;
                std::cout << "Enter your choice: ";
                std::getline(std::cin, input);
                trim(input);
                switch (input[0]) {
                    case '1': {
                        std::cout << "Enter the encrypted message: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class c(input);
                        mpz_class m;
                        mpz_powm(m.get_mpz_t(), c.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());
                        std::cout << "Decrypted message: " << m << std::endl;
                        crackLoop = false;
                        break;
                    }
                    case '2': {
                        std::cout << "Enter the encrypted message: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class c(input);
                        mpz_class m;
                        mpz_powm(m.get_mpz_t(), c.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());
                        std::cout << "Decrypted message: " << mpz_to_ascii_string(m) << std::endl;
                        crackLoop = false;
                        break;
                    }
                    default: {
                        std::cout << "Invalid choice" << std::endl;
                    }
                }
            }
            continue;
        }

        if (seq(input, "o") || seq(input, "other")) {
            bool otherLoop= true;
            while (otherLoop) {
                std::cout << "Choose an operation: " << std::endl;
                std::cout << "     [1] Big number calculator" << std::endl;
                std::cout << "     [2] String to int (ASCII)" << std::endl;
                std::cout << "     [3] int (ASCII) to String" << std::endl;
                std::cout << "     [4] String to binary" << std::endl;
                std::cout << "     [5] Binary to String" << std::endl;
                std::cout << "     [6] Convert between bases" << std::endl;
                std::cout << "     [7] Get big Primes" << std::endl;
                std::cout << "     [8] Get pair of Primes big enough for encrypting a specific message" << std::endl;
                std::cout << "     [9] get Big n for testing" << std::endl;
                std::cout << "     [D] find big dihedral Prime" << std::endl;
                std::cout << "     [P] check if number is Prime with high certainty" << std::endl;
                std::cout << "     [B] Back" << std::endl;
                std::cout << "     [Q] Quit" << std::endl;
                std::cout << "Enter your choice: ";
                std::getline(std::cin, input);
                trim(input);
                switch (input[0]) {
                    case '1': {
                        std::cout << "Enter the 1. Number: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class num1(input);
                        std::cout << "Enter the 2. Number: ";
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
                        break;
                    }
                    case '2': {
                        std::cout << "Enter a string: ";
                        std::getline(std::cin, input);
                        trim(input);
                        std::cout << "ASCII: " << ascii_string_to_mpz(input) << std::endl;
                        break;
                    }
                    case '3': {
                        std::cout << "Enter an int: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class num(input);
                        std::cout << "String: " << mpz_to_ascii_string(num) << std::endl;
                        break;
                    }

                    case '6': {
                        std::cout << "Useful Bases: " << std::endl;
                        std::cout << "     [2] Binary" << std::endl;
                        std::cout << "     [8] Octal" << std::endl;
                        std::cout << "     [10] Decimal" << std::endl;
                        std::cout << "     [16] Hexadecimal" << std::endl;
                        std::cout << "     [36] Base 36 (0-z)" << std::endl;
                        std::cout << "     [62] Base 62 (0-Z)" << std::endl;
                        std::cout << "Enter the base of the input(2-62): ";
                        std::getline(std::cin, input);
                        trim(input);
                        int from_base = std::stoi(input);
                        std::cout << "Enter the base of the output(2-62): ";
                        std::getline(std::cin, input);
                        trim(input);
                        int to_base = std::stoi(input);
                        std::cout << "Enter the number: ";
                        std::getline(std::cin, input);
                        trim(input);
                        std::cout << "Converted: " << convert_base(input, from_base, to_base) << std::endl;
                        break;
                    }
                    case '7': {
                        std::cout << "Enter the number from where to generate the primes: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class num(input);
                        std::cout << "The next 10 primes after " << num << " are:";
                        for (int i = 0; i < 10; i++) {
                            mpz_nextprime(num.get_mpz_t(), num.get_mpz_t());
                            std::cout << "\n     " << num << "\n";
                        }
                        std::cout << std::endl;
                        break;
                    }
                    case '8': {
                        std::cout << "Do you want to encrypt a string or an int?" << std::endl;
                        std::cout << "     [1] int" << std::endl;
                        std::cout << "     [2] string" << std::endl;
                        std::cout << "Enter your choice: ";
                        std::getline(std::cin, input);
                        trim(input);
                        switch (input[0]) {
                            case '1': {
                                std::cout << "Enter your message: ";
                                std::getline(std::cin, input);
                                trim(input);
                                mpz_class m(input);
                                gmp_randclass rng(gmp_randinit_default);
                                rng.seed(static_cast<unsigned long>(time(nullptr)));

                                mpz_class p;
                                mpz_class q;
                                mpz_class sqrt_m;
                                mpz_sqrt(sqrt_m.get_mpz_t(), m.get_mpz_t());

                                mpz_class delta;
                                mpz_class min_gap = sqrt_m / 20;
                                mpz_class max_gap = sqrt_m / 7;

                                while (true) {
                                    mpz_class p_guess = sqrt_m + 1 + rng.get_z_range(sqrt_m);
                                    mpz_nextprime(p.get_mpz_t(), p_guess.get_mpz_t());
                                    mpz_class q_min;
                                    mpz_cdiv_q(q_min.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());
                                    mpz_nextprime(q.get_mpz_t(), q_min.get_mpz_t());
                                    if (p != q) break;
                                }
                                std::cout << "prime p = " << p << std::endl;
                                std::cout << "prime q = " << q << std::endl;
                                std::cout << "Because " << p << " * " << q << " >= " << m << std::endl;
                            }
                            case '2': {
                                std::cout << "Enter your message: ";
                                std::getline(std::cin, input);
                                mpz_class m(ascii_string_to_mpz(input));

                                gmp_randclass rng(gmp_randinit_default);
                                rng.seed(static_cast<unsigned long>(time(nullptr)));

                                mpz_class p;
                                mpz_class q;
                                mpz_class sqrt_m;
                                mpz_sqrt(sqrt_m.get_mpz_t(), m.get_mpz_t());

                                mpz_class delta;
                                mpz_class min_gap = sqrt_m / 20;
                                mpz_class max_gap = sqrt_m / 7;

                                while (true) {
                                    mpz_class p_guess = sqrt_m + 1 + rng.get_z_range(sqrt_m);
                                    mpz_nextprime(p.get_mpz_t(), p_guess.get_mpz_t());
                                    mpz_class q_min;
                                    mpz_cdiv_q(q_min.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());
                                    mpz_nextprime(q.get_mpz_t(), q_min.get_mpz_t());
                                    if (p != q) break;
                                }
                                std::cout << "prime p = " << p << std::endl;
                                std::cout << "prime q = " << q << std::endl;
                                std::cout << "Because " << p << " * " << q << " >= " << m << " = " << input << std::endl;
                            }
                            default: std::cout << "Enter a valid choice" << std::endl;
                        }
                        break;
                    }
                    case '9': {
                        std::cout << "How big should n be?: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class n(input);
                        mpz_sqrt(n.get_mpz_t(), n.get_mpz_t());
                        mpz_nextprime(n.get_mpz_t(), n.get_mpz_t());
                        mpz_class m(n);
                        mpz_nextprime(n.get_mpz_t(), n.get_mpz_t());
                        std::cout << "n = " << n*m << std::endl;
                        break;
                    }
                    case 'd':
                    case 'D': {
                        std::cout << "how big should the prime be at least?: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class n(input);

                        mpz_class candidate(n);

                        while (true) {
                            candidate = next_palindrome(candidate);
                            if(mpz_probab_prime_p(candidate.get_mpz_t(), 100)) {
                                std::cout << "Found prime: " << candidate << std::endl;
                                break;
                            }
                        }
                        break;
                    }
                    case 'p':
                    case 'P': {
                        std::cout << "Enter Number: ";
                        std::getline(std::cin, input);
                        trim(input);
                        mpz_class n(input);
                        if (mpz_probab_prime_p(n.get_mpz_t(), 1000) > 0) {
                            std::cout << "Number is prime" << std::endl;
                        }
                        else {
                            std::cout << "Number is not prime" << std::endl;
                        }
                        break;
                    }
                    case 'b':
                    case 'B':
                        otherLoop = false;
                    break;
                    case 'q':
                    case 'Q':
                        mainloop = false;
                    otherLoop = false;
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

