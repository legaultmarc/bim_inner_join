#include <unordered_set>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>


struct Variant {
    unsigned int chrom;
    std::string name;
    unsigned int position;
    std::string a1;
    std::string a2;

    bool operator==(const Variant& that) const {
        return chrom == that.chrom && position == that.position;
    };

    bool operator!=(const Variant& that) const {
        return !(*this == that);
    };

    bool operator<(const Variant& that) const {
        if (chrom != that.chrom)
            return chrom < that.chrom;
        else
            return position < that.position;
    };

    bool operator>(const Variant& that) const {
        if (chrom != that.chrom)
            return chrom > that.chrom;
        else
            return position > that.position;
    };

    bool alleles_eq(const Variant& that) const {
        bool locus_match = *this==that;

        // Allele comparison is made a bit harder because sometimes there are
        // "0" in bim files if the variant is monomorphic.
        // This means that A/0 and G/A should match.
        // On the other hand, something like A/0 and T/G should not.
        // One way of testing this would be to check that the number of unique
        // non-"0" alleles is <= 2.

        std::unordered_set<std::string> alleles;
        if (a1 != "0")
            alleles.insert(a1);
        if (a2 != "0")
            alleles.insert(a2);
        if (that.a1 != "0")
            alleles.insert(that.a1);
        if (that.a2 != "0")
            alleles.insert(that.a2);

        bool alleles_match = alleles.size() <= 2;

        return locus_match && alleles_match;
    };
};


struct Log {
    int n_files;
    std::ofstream* names_lists;
    std::ofstream matches_bim;
    std::ofstream mismatches_bim;

    void open(int n_files) {
        _open = true;

        // Open the output files for the lists.
        this->n_files = n_files;
        names_lists = new std::ofstream[n_files];

        for (int i = 0; i < n_files; i++) {
            std::ostringstream filename;
            filename << "bij_names_" << i + 1 << ".txt";
            names_lists[i].open(filename.str());

            if (!names_lists[i].is_open()) {
                std::cerr << "Could not write to "
                          << filename.str() << std::endl;
                throw 1;
            }
        }

        // Open the output file for the matches.
        matches_bim.open("bij_matches.bim");
        if (!matches_bim.is_open()) {
            std::cerr << "Could not write to " << "bij_matches.bim" << std::endl;
            throw 1;
        }

        mismatches_bim.open("bij_mismatches.bim");
        if (!mismatches_bim.is_open()) {
            std::cerr << "Could not write to "
                      << "bij_mismatches.bim" << std::endl;
            throw 1;
        }
    };

    ~Log() {
        if (_open) {
            for (int i = 0; i < n_files; i++) {
                names_lists[i].close();
            }
            delete[] names_lists;

            matches_bim.close();
            mismatches_bim.close();
        }
    };

    private:
        bool _open = false;

};



// Output string formatting for variants.
std::ostream& operator<<(std::ostream &s, const Variant &v) {
    return s
        << "<Variant " << v.name
        << " chr" << v.chrom << ":" << v.position
        << ", [" << v.a1 << ", " << v.a2 << "]" << ">";
}


int print_usage() {
    std::cout << "Usage: " << std::endl;
    std::cout << "\tbim-inner-join file1.bim [ file2.bim, ... ]" << std::endl;
    return 1;
}


// Reads a line from a bim file and returns a Variant instance.
bool read_variant(std::ifstream &file, Variant &v) {
    std::string line;
    getline(file, line);

    if (line == "")
        return false;

    std::istringstream ss(line);
    std::string cm;

    ss >> v.chrom >> std::ws
       >> v.name >> std::ws
       >> cm >> std::ws
       >> v.position >> std::ws
       >> v.a1 >> std::ws
       >> v.a2;

    return true;
}


void print_variants(const Variant* li, int n) {
    for (int i = 0; i < n - 1; i++) {
        std::cout << li[i] << ", ";
    }
    std::cout << li[n - 1] << std::endl;
}


const Variant* max(const Variant* li, int n) {
    const Variant *greatest = li;
    for (int i = 1; i < n; i++) {
        if (li[i] > *greatest) {
            greatest = li + i;
        }
    }
    return greatest;

}

// Read variants from the files if they are before the maximum.
void read_variants_until_max(
    std::ifstream *files,
    Variant *cur,
    const Variant& max,
    const int n
    ) {

    for (int i = 0; i < n; i++) {
        if (cur[i] < max && !files[i].eof()) {
            read_variant(files[i], cur[i]);
        }
    }
}


void read_variants(std::ifstream *files, Variant *cur, int n) {
    for (int i = 0; i < n; i++) {
        if (!files[i].eof())
            read_variant(files[i], cur[i]);
    }
}


bool any_eof(const std::ifstream *files, int n) {
    for (int i = 0; i < n; i++) {
        if (files[i].eof())
            return true;
    }
    return false;
}


void write_bim_line(const Variant &v, std::ofstream &file) {
    file << v.chrom << "\t"
         << v.name << "\t"
         << "0" << "\t"
         << v.position << "\t"
         << v.a1 << "\t"
         << v.a2 << std::endl;
}


// Check if all the variants in the list are equals.
bool check_match(const Variant* li, int n, Log &log) {
    bool wrote_bim_line = false;
    bool all_match = true;

    const Variant first = li[0];
    for (int i = 1; i < n; i++) {
        if (!li[i].alleles_eq(first))
            all_match = false;
    }

    if (all_match) {
        // Write the names to the appropriate lists and to the "matches" bim
        // file.
        for (int i = 0; i < n; i++) {
            log.names_lists[i] << li[i].name << std::endl;

            // Find a bim line with no "0" allele to write, otherwise
            // write the first one.
            if (!wrote_bim_line && li[i].a1 != "0" && li[i].a2 != "0") {
                write_bim_line(li[i], log.matches_bim);
                wrote_bim_line = true;
            }
        }

        if (!wrote_bim_line)
            write_bim_line(li[0], log.matches_bim);
    }
    else {
        // TODO I am not sure how to log the mismatches...
        // write_bim_line(first, log.mismatches_bim);
    }
    // print_variants(li, n);

    return all_match;
}


int main(int argc, char* argv[]) {
    int n_files = argc - 1;
    std::ifstream files[n_files];

    Log log;
    try {
        log.open(n_files);
    }
    catch (int code) {
        return 1;
    }

    // Open file pointers to the files to join and check that they exist and
    // that at least two files were provided.
    if (argc <= 2) {
        return print_usage();
    }

    for (int i = 1; i < argc; i++) {
        std::cout << "Opening: " << argv[i] << std::endl;
        files[i - 1].open(argv[i]);

        if (!files[i - 1].is_open()) {
            std::cout << "Could not find file: " << argv[i] << std::endl;
            return 1;
        }
    }

    Variant cur[n_files];
    read_variants(files, cur, n_files);

    while (!any_eof(files, n_files)) {
        if (!check_match(cur, n_files, log)) {
            // Find the furthest variant.
            const Variant *max_variant = max(cur, n_files);
            // std::cout << "Maximum: " << *max_variant << std::endl;

            // Step the others towards the maximum.
            read_variants_until_max(files, cur, *max_variant, n_files);

            // std::cout << "Now : ";
            // print_variants(cur, n_files);
        }

        else {
            read_variants(files, cur, n_files);
        }
    }

    return 0;
}
