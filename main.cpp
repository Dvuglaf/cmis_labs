#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <cmath>
#include "fftw3.h"
#include "boost/math/special_functions/gamma.hpp"
#include <fstream>

const std::vector<bool> c_polynom_first  = { 1, 1, 1, 1, 0, 0, 0, 1 };
const std::vector<bool> c_polynom_second = { 1, 1, 0, 0, 0, 0, 0, 1 };
const std::vector<bool> c_polynom_third  = { 1, 0, 0, 1, 0, 0, 0, 1 };


// Период последовательности [begin; end)
template <typename RandomAccessIt>
size_t period(RandomAccessIt begin, RandomAccessIt end) {
	const size_t seq_size = std::distance(begin, end);
	
	// перебор возможного периода от 1 до размера последовательности, возвращается минимально подходящий период
	for (size_t subseq_size = 1; subseq_size < seq_size; ++subseq_size) {
		// проверка, является ли subseq_size периодом
		for (size_t j = subseq_size; j < seq_size; ++j) {
			if (*(begin + j) == *(begin + (j % subseq_size))) {
				if (j == seq_size - 1) return subseq_size; // если дошли до последней итераци
			}
			else break;
		}
	}
	return seq_size;
}


// Один такт LFSR; a - коэффициенты полинома, S - состояния
bool lfsr(const std::vector<bool>& a, std::vector<bool>& S) {
	const auto N = a.size();
		
	auto acc = false;		// накапливаемая сумма
	const bool ret = S[0];	// первый элемент, его возвращаем (до свдига)
	
	// вычисление суммы
	for (size_t i = 0; i < N; ++i) {
		if (a[i] == 1) {
			acc ^= S[i];
		}
	}

	// сдвиг
	for (size_t i = 0; i < N - 1; ++i) {
		S[i] = S[i + 1];
	}
	S[N - 1] = acc;

	return ret;
}

// Комбинация двух LFSR, length - требуемая длина генерируемой последовательности
std::vector<bool> sg_generator(const size_t length) {
	std::vector<bool> sequence;
	sequence.reserve(length);

	// Полиномы
	const std::vector<bool> polynom_first = { 1, 1, 0, 0, 0, 0, 0 };
	const std::vector<bool> polynom_second = { 1, 1, 0, 0, 0, 0 };
	
	// Начальные состояния
	std::vector<bool> S_first = { 1, 0, 1, 0, 0, 1, 1 };
	std::vector<bool> S_second = { 1, 1, 0, 1, 0, 1 };

	// Генерация последовательности
	while (sequence.size() < length) {
		auto a_i = lfsr(polynom_first, S_first);
		auto s_i = lfsr(polynom_second, S_second);
		// Проверка на бит '1' селектирующей последовательности
		if (s_i) sequence.push_back(a_i);
	}
	return sequence;
}


bool spectral_test(const std::vector<bool>& eps, const size_t N) {
	fftw_complex* in = nullptr, *out = nullptr;
	fftw_plan my_plan = nullptr;
	double* seq_arr = nullptr;
	try {
		in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
		out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
		seq_arr = new double[N];
		for (size_t i = 0; i < N; ++i) {
			seq_arr[i] = eps[i] ? 1 : -1;
		}
		my_plan = fftw_plan_dft_r2c_1d(N, seq_arr, out, FFTW_ESTIMATE);
		fftw_execute(my_plan);
		// Массив модулей комплексных чисел
		std::vector<double> mod;
		mod.reserve(N / 2 + 1);
		for (size_t i = 0; i < std::ceil(N / 2); ++i) {
			mod.push_back((std::sqrt((out[i][0] * out[i][0]) + (out[i][1] * out[i][1]))));
		}
		auto T = std::sqrt(N * std::log(1 / 0.05));
		auto N_0 = 0.95 * N / 2;
		auto N_1 = 0;
		for (auto m : mod) {
			if (m < T) ++N_1;
		}
		auto d = (N_1 - N_0) / std::sqrt((N * (0.95 * 0.05) / 4));
		auto P_value = std::erfc(std::fabs(d) / std::sqrt(2));
		std::cout << "P_value = " << P_value << std::endl;

		fftw_destroy_plan(my_plan);
		fftw_free(in);
		fftw_free(out);
		delete[] seq_arr;

		if (P_value > 0.01) return true;
		return false;
	}
	catch (std::exception& e) {
		fftw_destroy_plan(my_plan);
		fftw_free(in);
		fftw_free(out);
		delete[] seq_arr;
		std::cout << e.what();
		return false;
	}
}

size_t slider(const std::vector<bool>& block, const std::vector<bool>& B) {
	size_t count = 0;
	size_t k = 0; // Число совппадение окна и слайдера
	for (size_t i = 0; i <= block.size() - B.size(); ++i) { // Такое условие чтобы не выйти за границы
		for (size_t j = 0; j < B.size(); ++j) {
			if ((k == B.size() - 1) && (block[i + j] == B[j])) {
				// Здесь только если все биты окна и слайдера совпали
				++count;
				k = 0;
			}
			else if (block[i + j] == B[j]) {
				++k;
			}  // В этом случае дальнейшая проверка следующего бита
			else {
				k = 0;
				break;
			}
		}
	}
	return count;
}

auto partition(const std::vector<bool>& eps, const size_t M, const size_t N) {
	std::vector<std::vector<bool>> ret;
	ret.reserve(M);
	size_t k = 0; // Индекс для обхода eps
	for (size_t i = 0; i < N; ++i) { // Создадим N блоков
		std::vector<bool> block;
		block.reserve(N);
		for (size_t j = 0; j < M; ++j) { // Переместим M элементов в блок
			block.push_back(eps[k]);
			++k;
		}
		ret.push_back(block);
	}
	return ret;
}

bool overlapping_templates_test(std::vector<bool>& eps, const std::vector<bool>& B,
						   const size_t M = 1000, const size_t N = 1000) {
	auto test = partition(eps, M, N);

	std::vector<size_t> v = { 0, 0, 0, 0, 0, 0 };
	for (auto& block : test) {
		auto cnt = slider(block, B);
		if (cnt < 6) {
			v[cnt] += 1;
		}
		if (cnt >= 6)
			v[5] += 1;
	}
	double xi_s = 0.0;
	const std::vector<double> pi = { 0.364091, 0.185659, 0.139381,
									 0.100571, 0.070432, 0.139865 };
	
	for (size_t i = 0; i < pi.size(); ++i) {
		xi_s += (std::pow(v[i] - N * pi[i], 2) / (N * pi[i]));
	}
	std::cout << "v_i = ";
	for (auto& el : v) {
		std::cout << el << ' ';
	}
	std::cout << "\nxi_s = " << xi_s << std::endl;
	auto P_value = boost::math::gamma_q<double>(2.5, xi_s/2);
	std::cout << "P_value = " << P_value << std::endl;

	if (P_value > 0.01) return true;
	return false;
}


size_t linear_complexity(const std::vector<bool>& s) {
	const auto n = s.size();
	std::vector<bool> b(n, false);
	std::vector<bool> t(n, false);
	std::vector<bool> c(n, false);

	b[0] = true;
	c[0] = true;
	int N = 0;
	int L = 0;
	int m = -1;

	while (N < n) {
		bool d = s[N];
		for (size_t i = 1; i <= L; ++i) {
			d ^= c[i] & s[N - i];
		}
		if (d) {
			//t = c;
			for (size_t i = 0; i < n; ++i) {
				t[i] = c[i];
			}
			for (size_t i = 0; (i + N - m) < n; ++i) {
				c[i + N - m] = c[i + N - m] ^ b[i];
			}
			if (2 * L <= N) {
				L = N + 1 - L;
				m = N;
				//b = t;
				for (size_t i = 0; i < n; ++i) {
					b[i] = t[i];
				}
			}
		}
		++N;
	}
	return L;

}


bool linear_complexity_test(const std::vector<bool>& eps, const size_t M) {
	const size_t N = eps.size() / M;
	auto mu = M / 2. + (9 + std::pow(-1, M + 1)) / 36 + (M / 3. + 2 / 9) / std::pow(2, M);
	auto blocks = partition(eps, M, N);
	std::vector<double> T;
	T.reserve(N);
	for (auto& block : blocks) {
		auto L_i = linear_complexity(block);
		T.push_back(std::pow(-1, M) * (L_i - mu) + 2 / 9);
	}
	std::vector<size_t> v(7, 0);
	for (auto& T_i : T) {
		if (T_i <= -2.5) ++v[0];
		else if (T_i > -2.5 && T_i <= -1.5) ++v[1];
		else if (T_i > -1.5 && T_i <= -0.5) ++v[2];
		else if (T_i > -0.5 && T_i <= 0.5) ++v[3];
		else if (T_i > 0.5 && T_i <= 1.5) ++v[4];
		else if (T_i > 1.5 && T_i <= 2.5) ++v[5];
		else if (T_i > 2.5) ++v[6];
	}
	std::cout << "v_i = ";
	for (auto v_i : v) {
		std::cout << v_i << ' ';
	}
	std::cout << std::endl;

	const std::vector<double> pi = { 0.010417, 0.03125, 0.125, 0.5,
									 0.25, 0.0625, 0.020833 };
	double xi_s = 0.0;
	for (size_t i = 0; i < pi.size(); ++i) {
		xi_s += (std::pow(v[i] - N * pi[i], 2) / (N * pi[i]));
	}
	std::cout << "xi_s = " << xi_s << std::endl;
	auto P_value = boost::math::gamma_q<double>(7 / 2, xi_s / 2);
	std::cout << "P_value = " << P_value << std::endl;

	if (P_value > 0.01) return true;
	return false;


}

int main() {
	constexpr size_t N = 1000000;
	std::vector<bool> seq = sg_generator(N);
	std::cout << "Period of sequence = " << period(seq.begin(), seq.end()) << std::endl;

	// tests for e
	{
		std::fstream file("C:\\Users\\aizee\\Downloads\\2\\data.e", std::fstream::in);
		std::vector<bool> e;
		e.reserve(N);

		if (file.is_open()) {
			char value = 0;
			while (!file.eof()) {
				file >> value;
				e.push_back(value == '1' ? true : false);
				if (e.size() == N)
					break;
			}

		}
		file.close();
		std::cout << "Spectral test e: \n";
		spectral_test(e, e.size());

		std::cout << "\nOverlapping test e: \n";
		overlapping_templates_test(e, { 1, 1, 1, 1, 1, 1, 1, 1, 1 }, 1032, 968);

		std::cout << "\nLinear test e: \n";
		linear_complexity_test(e, 500);
	}

	//tests for pi
	{
		std::fstream file("C:\\Users\\aizee\\Downloads\\2\\data.pi", std::fstream::in);
		std::vector<bool> pi;
		pi.reserve(N);

		if (file.is_open()) {
			char value = 0;
			while (!file.eof()) {
				file >> value;
				pi.push_back(value == '1' ? true : false);
				if (pi.size() == N)
					break;
			}

		}
		file.close();
		std::cout << "\nSpectral test pi : \n";
		spectral_test(pi, pi.size());

		std::cout << "\nOverlapping test pi: \n";
		overlapping_templates_test(pi, { 1, 1, 1, 1, 1, 1, 1, 1, 1 }, 1032, 968);

		std::cout << "\nLinear test pi: \n";
		linear_complexity_test(pi, 500);
	}

	std::cout << "\n______________________________________\n";
	std::cout << "Spectral test:\n";
	if (spectral_test(seq, N)) std::cout << "Sequence is random\n";
	else std::cout << "Sequence is NOT random\n";
	std::cout << std::endl;

	std::cout << "Overlapping templates test:\n";
	if (overlapping_templates_test(seq, { 1, 1, 1, 1, 1, 1, 1, 1, 1 })) std::cout << "Sequence is random\n";
	else std::cout << "Sequence is NOT random\n";
	std::cout << std::endl;

	std::cout << "Linear complexity test::\n";
	if (linear_complexity_test(seq, 500)) std::cout << "Sequence is random\n";
	else std::cout << "Sequence is NOT random\n";
	std::cout << std::endl;
	
}


//{ 0.324652, 0.182617, 0.142670, 0.106645, 0.077147, 0.166269 };
