#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <map>
#include <cstdio>
using namespace std;

vector<string> str_split(string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != string::npos) {
        token = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
}

int int_from(string s) {
    if (s[0] == '-') {
        return -1 * int_from(s.substr(1));
    }

    int n = 0;
    for (int i = 0; i < s.length(); i++) {
        n *= 10;
        n += (s[i] - '0');
    }

    return n;
}

map<pair<string, string>, int> counter;
map<pair<string, string>, string> sequence_map;
map<pair<string, string>, vector<int> > tensor_map;
map<pair<string, string>, vector<double> > probabilities_map;

map<pair<string, string>, int>::iterator it;
map<pair<string, string>, int>::iterator counter_iterator;
map<pair<string, string>, string>::iterator sequence_map_iterator;
map<pair<string, string>, vector<int> >::iterator tensor_map_iterator;
map<pair<string, string>, vector<double> >::iterator probabilities_map_iterator;

int main(int argc, char* argv[]) {
    int minimum_count_to_output = atoi(argv[1]);
    string s;
    int i;

    while (getline(cin, s)) {
        vector<string> v = str_split(s, "\t");
        string chromosome = v[0];
        string position = v[1];
        string sequence = v[2];

        pair<string, string> key = make_pair(chromosome, position);

        counter[key] += 1;

        sequence_map_iterator = sequence_map.find(key);
        if (sequence_map_iterator == sequence_map.end()) {
            sequence_map[key] = sequence;
        }

        tensor_map_iterator = tensor_map.find(key);
        if (tensor_map_iterator == tensor_map.end()) {
            vector<string> tensor_str(v.begin()+3, v.begin()+3 + 33*8*4);

            vector<int> tensor;
            for (i = 0; i < tensor_str.size(); i++) {
                tensor.push_back(int_from(tensor_str[i]));
            }
            tensor_map[key] = tensor;
        }

        vector<string> probabilities_str(v.begin()+3 + 33*8*4, v.end());
        vector<double> probabilities;
        for (i = 0; i < probabilities_str.size(); i++) {
            probabilities.push_back(atof(probabilities_str[i].c_str()));
        }

        probabilities_map_iterator = probabilities_map.find(key);
        if (probabilities_map_iterator != probabilities_map.end()) {
            vector<double> p = probabilities_map_iterator->second;
            for (i = 0; i < probabilities.size(); i++) {
                p[i] = p[i] + probabilities[i];
            }
            probabilities_map[key] = p;
        } else {
            probabilities_map[key] = probabilities;
        }
    }

    for (it = counter.begin(); it != counter.end(); it++) {
        pair<string, string> key = it->first;
        int count = it->second;
        if (count < minimum_count_to_output) {
            continue;
        }

        string chromosome = key.first;
        string position = key.second;
        string sequence = sequence_map[key];
        vector<int> tensor = tensor_map[key];
        vector<double> probabilities = probabilities_map[key];

        printf("%s\t%s\t%s", chromosome.c_str(), position.c_str(), sequence.c_str());
        for (i = 0; i < tensor.size(); i++) {
            printf("\t%d", tensor[i]);
        }
        for (i = 0; i < probabilities.size(); i++) {
            printf("\t%.6f", probabilities[i] / (count + 0.0));
        }

        printf("\n");
    }

    return 0;
}
