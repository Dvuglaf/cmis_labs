#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <random>
#include <fstream>
#include <string>
#include <regex>
#include <map>

std::vector<std::string> permutations;

void generate_permutations(const std::string& alphabet, const std::string& prefix, const int N, const int n) {
    // base case: prefix + all symbols of alphabet
    if (n == 1) {
        for (int i = 0; i < N; ++i)
            permutations.push_back(prefix + alphabet[i]);    
    }
    // recursion case
    else {
        for (int i = 0; i < N; ++i)
            generate_permutations(alphabet, prefix + alphabet[i], N, n - 1);
    }
}

void save_table_to_file(const std::unordered_map<std::string, std::string>& table, const std::string& name_without_extention) {
    // path = name_without_extention_n.txt, where n is number of dimension
    const std::string path = name_without_extention + "_" + static_cast<char>(table.begin()->first.size() + 48) + ".txt";

    std::ofstream out; 
    out.open(path);
    if (out.is_open()) {
        for (const auto& it : table) {

            // write to file pair KEY:VALUE on each line 
            out << it.first << ":" << it.second << "\n";

        }
    }
    out.close();
}

// split string by delimiter, result is array of substrings
std::vector<std::string> split(const std::string& string, const std::string& delimiter) {
    std::vector<std::string> tokens;
    size_t previous = 0, position = 0;
    do {
        position = string.find(delimiter, previous);
        if (position == std::string::npos) position = string.length();
        std::string token = string.substr(previous, position - previous);
        if (!token.empty()) tokens.push_back(token);
        previous = position + delimiter.length();
    } while (position < string.length() && previous < string.length());
    return tokens;
}

// split string into blocks of fix size = n
std::vector<std::string> split(const std::string& string, const int n) {
    std::vector<std::string> tokens;
    size_t previous = 0, position = 0;
    do {
        position += n;
        if (position == std::string::npos) position = string.length();
        std::string token = string.substr(previous, position - previous);
        if (!token.empty()) tokens.push_back(token);
        previous = position;
    } while (position < string.length() && previous < string.length());
    return tokens;
}

std::unordered_map<std::string, std::string> read_table_from_file(const std::string& path) {
    std::unordered_map<std::string, std::string> table;
    std::ifstream in;
    std::string line;
    in.open(path);
    if (in.is_open()) {
        while (std::getline(in, line)) {
            // splitted = [key, value]
            const auto splitted = split(line, ":");
            table.insert(std::make_pair(splitted[0], splitted[1]));
        }
    }
    in.close();
    return table;
}

std::string read_strings_from_file(const std::string& path) {
    std::ifstream in;
    std::string line, string;
    in.open(path);
    if (in.is_open()) {
        while (std::getline(in, line)) {
            string += line;
            string += '\n';
        }
        // delete last extra '\n' symbol
        string.erase(--string.end());
    }
    in.close();
    return string;
}

void save_string_to_file(const std::string& string, const std::string& path) {
    std::ofstream out;      
    out.open(path); 
    if (out.is_open()) {
        out << string;
    }
    out.close();
}

std::unordered_map<std::string, std::string> generate_table(const std::string& alphabet, const int N, const int n) {
    generate_permutations(alphabet, "", N, n);

    std::random_device rd;
    std::mt19937 gen(rd());

    // array of random indexes, will be used to fill the map
    std::vector<int> random_shuffled;
    random_shuffled.reserve(permutations.size());
    for (int i = 0; i < permutations.size(); ++i) {
        random_shuffled.push_back(i);
    }
    std::shuffle(random_shuffled.begin(), random_shuffled.end(), gen);

    std::unordered_map<std::string, std::string> table;
    table.reserve(permutations.size());

    int i = 0;
    for (const auto& it : permutations) {
        table[it] = permutations[random_shuffled[i]];
        //table.insert(std::make_pair(it, permutations[random_shuffled[i]]));
        ++i;
    }
    permutations.clear();
    return table;
}

void save_decipher_table(const std::unordered_map<std::string, std::string>& cipher_table, const std::string& name_without_extention) {
    // path = name_without_extention_n.txt, where n is number of dimension
    const std::string path = name_without_extention + "_" + static_cast<char>(cipher_table.begin()->first.size() + 48) + ".txt";

    std::ofstream out;
    out.open(path);
    if (out.is_open()) {
        for (const auto& it : cipher_table) {

            // write to file pair VALUE:KEY on each line 
            out << it.second << ":" << it.first << "\n";

        }
    }
    out.close();
}


// save indexes of symbols ' ' and '\n' from input string
std::map<int, char> index_of_deleted_symbols(const std::string& string) {
    std::map<int, char> indexes;
    size_t prev = 0;
    while (string.find(' ', prev) != std::string::npos) {
        const auto pos = string.find(' ', prev);
        indexes.insert(std::make_pair(pos, ' '));
        prev = pos + 1;
    }
    prev = 0;
    while (string.find('\n', prev) != std::string::npos) {
        const auto pos = string.find('\n', prev);
        indexes.insert(std::make_pair(pos, '\n'));
        prev = pos + 1;
    }
    return indexes;
}

std::string cipher_or_decipher(const std::unordered_map<std::string, std::string> table, const std::string& string, const int n) {
    const auto indexes_of_deleted_symbols = index_of_deleted_symbols(string);
    std::string replaced_string = std::regex_replace(string, std::regex("\[ ,\n]"), "");
    const auto splitted = split(replaced_string, n);
    std::string result;
    
    for (const auto& it : splitted) {
        if (it.size() != n) {
            result += it;
            continue;
        }
        result += table.at(it);
    }
    for (const auto& it : indexes_of_deleted_symbols) {
        result.insert(it.first, 1, it.second);
    }
    
    return result;
}

int main() {
    
    // n- dimensions
    const int n = 4;

    // if generate equals true the table of permutations will be generated again
    const bool generate = false;

    // define cipher or decipher programm wil be used
    const bool cipher  = false;

    // symbols of english alphabet
    const std::string alphabet = "abcdefghijklmnopqrstuvwxyz";
    const auto N = alphabet.size();
    
    if (generate) {
        const auto table = generate_table(alphabet, N, n);
        save_table_to_file(table, "cipher_table");
        save_decipher_table(table, "decipher_table");
    }

    if (cipher) {
        const auto cipher_table = read_table_from_file("cipher_table_4.txt");

        auto input = read_strings_from_file("input.txt");
        
        std::cout << "Input string:  " << input << std::endl;
        const auto result = cipher_or_decipher(cipher_table, input, n);
        std::cout << "Result string: " << result << std::endl;
        save_string_to_file(result, "result.txt");

    }
    else {
        const auto decipher_table = read_table_from_file("decipher_table_4.txt");

        const auto input = read_strings_from_file("result.txt");
        std::cout << "Input string:  " << input << std::endl;
        const auto result = cipher_or_decipher(decipher_table, input, n);
        std::cout << "Result string: " << result << std::endl;
        save_string_to_file(result, "check.txt");

    }
    return 0;

}

