#include "gmp.h"
#include <iostream>
#include <vector>
#include "mpirxx.h"
#include "poly.h"

#include <NTL/ZZ_pX.h>
using namespace NTL;




const std::vector<int> PRIMES = {	2,3,5,7,11,13,17,19,23,29,31,37,
									41,43,47,53,59,61,67,71,73,79,83,89,
									97,101,103,107,109,113,127,131,137,139,149,151,
									157,163,167,173,179,181,191,193,197,199,211,223,
									227,229,233,239,241,251,257,263,269,271,277,281,
									283,293,307,311,313,317,331,337,347,349,353,359,
									367,373,379,383,389,397,401,409,419,421,431,433,
									439,443,449,457,461,463,467,479,487,491,499,503,
									509,521,523,541,547,557,563,569,571,577,587,593,
									599,601,607,613,617,619,631,641,643,647,653,659,
									661,673,677,683,691,701,709,719,727,733,739,743,
									751,757,761,769,773,787,797,809,811,821,823,827,
									829,839,853,857,859,863,877,881,883,887,907,911,
									919,929,937,941,947,953,967,971,977,983,991,997
								};

bool is_prime(const mpz_class& n) {
	if (n < 2) return false;
	if (n < 4) return true;
	if (n % 2 == 0) return false;
	mpz_class root;
	mpz_sqrt(root.get_mpz_t(), n.get_mpz_t());
	root += 1;
	for (mpz_class i = 3; i < root; i+=2) {
		if (n % i == 0) return false;
	}
	return true;
}

mpz_class max_prime_factor(const mpz_class& n) {
	mpz_class r = n;
	mpz_class max_fact = 0;
	if (r % 2 == 0) {
		max_fact = 2;
		do {
			r /= 2;
		} while (r == 2);
	}
	for (mpz_class i = 3; i <= r; i += 2) {
		if (r % i == 0) {
			max_fact = i;
			do {
				r /= i;
			} while (r == i);
		}
	}
	return max_fact;
}

bool aks(const mpz_class& n) {
	mpz_class q, r;
	// Check if n is perfect prime
	for (auto& prime : PRIMES) {
		mpz_fdiv_qr_ui(q.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t(), prime);
		if (r == 0 && q != 1) return false;
		if (prime > n) return true;
	}

	// Find minimal order r such ord_r(n) > (log(n))^2
	const mpz_class logn = mpz_sizeinbase(n.get_mpz_t(), 2);
	r = 2;
	mpz_class gcd;
	while (r < n) {
		mpz_gcd(gcd.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t());
		if (gcd != 1) return false;
		if (is_prime(r)) {
			q = max_prime_factor(r - 1);
			mpz_class temp;
			mpz_sqrt(temp.get_mpz_t(), r.get_mpz_t()); // temp = sqrt(r)
			if (q > 4 * temp * logn) {
				temp = (r - 1) / q;
				mpz_powm(temp.get_mpz_t(), n.get_mpz_t(), temp.get_mpz_t(), r.get_mpz_t());
				if (temp != 1) break;
			}
		}
		++r;
	}

	// 3.
	mpz_class min = (r < n - 1) ? r : n - 1;
	for (mpz_class a = 2; a <= min; ++a) {
		if (n % a == 0) return false;
	}

	// 4.
	if (n <= r) return true;

	// 5.
	mpz_class root_r;
	mpz_sqrt(root_r.get_mpz_t(), r.get_mpz_t());
	mpz_class check = 2 * root_r * logn;
	if (n - 1 <= check) {
		for (mpz_class a = r; a <= n - 1; ++a) {
			mpz_gcd(gcd.get_mpz_t(), a.get_mpz_t(), n.get_mpz_t());
			if (gcd != 1) return false;
		}
	}
	else {
		ZZ_pPush push(to_ZZ(n.get_str().c_str()));
		ZZ_p::init(to_ZZ(n.get_str().c_str())); // initialize GF(n)


		ZZ_pX left_poly, right_poly; // it will be (x+a)^n and x^n
		ZZ_pXModulus Mod_poly; // pre-computed information about x^r-1 to speed up computations, built from ZZ_pX
		ZZ_pX mod_poly; // it will be x^r-1

		mod_poly = ZZ_pX(r.get_ui(), 1); // 1 * x^r
		sub(mod_poly, mod_poly, 1); // x^r-1
		build(Mod_poly, mod_poly); // pre-compute information about x^r-1 and speed up
		
		std::cout << "\nNum of iterations = " << check.get_str() << std::endl;

		ZZ_p z_check = to_ZZ_p(check.get_ui());
		ZZ z_n = to_ZZ(n.get_str().c_str());
		std::cout << "r = " << r.get_str() << std::endl;
		// all computation over GF(n)
		for (ZZ_p a = to_ZZ_p(1); a != z_check; ++a) {
			PowerXPlusAMod(left_poly, a, z_n, Mod_poly); // left_poly = (x+a)^n % x^r-1
			PowerXMod(right_poly, z_n, Mod_poly); //right_poly = x^n % x^r-1
			
			if (!(IsZero((left_poly - right_poly - a) % Mod_poly))) return false; //((x+a)^n - x^n - a) % x^r-1
			std::cout << a << std::endl;

		}
	}
	return true;

}


int main() {
	mpz_class a(std::string("2147483647"));
	mpz_class b(std::string("62813549"));
	std::cout << std::boolalpha << aks(a * b);
}
// mpz_class mersenn(std::string("2305843009213693951"));
