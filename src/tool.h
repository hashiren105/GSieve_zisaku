#ifndef TOOL_H
#define TOOL_H

#pragma once
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stack>
#include <list>
#include <cctype>

using namespace std;
using namespace NTL;

bool read_basis(const string& input_path, mat_ZZ& B);
bool write_basis(const mat_ZZ& M, const string& filename);

vec_RR to_RR(const vec_ZZ& v);
mat_RR to_RR(const mat_ZZ& M);

RR computeGH(const mat_ZZ &B);
RR Computenorm(const vec_ZZ &v);

bool Reduce(vec_ZZ* p1, const vec_ZZ& p2);

// Template function declarations

template <typename T>
bool read_basis(const string& input_path, Mat<T>& B) {
    ifstream in(input_path);
    if (!in) {
        cerr << "Failed to open input file: " << input_path << endl;
        return false;
    }

    string line;
    vector<vector<T>> temp_rows;

    while (getline(in, line)) {
        vector<T> row;
        string num;
        string trimmed = line;
        size_t start = trimmed.find_first_not_of(" \t\r\n");
        if (start == string::npos)
            continue;
        size_t end = trimmed.find_last_not_of(" \t\r\n");
        trimmed = trimmed.substr(start, end - start + 1);

        if (trimmed == "]")
            break;
        for (char ch : line) {
            if (isdigit(ch) || ch == '-' || ch == '+') {
                num += ch;
            } else if (!num.empty()) {
                row.push_back(conv<T>(num.c_str()));
                num.clear();
            }
        }
        if (!num.empty()) {
            row.push_back(conv<T>(num.c_str()));
        }
        if (!row.empty()) {
            if (!temp_rows.empty() && row.size() != temp_rows[0].size())
                continue;
            temp_rows.push_back(row);
        }
    }

    if (temp_rows.empty()) {
        cerr << "Input file " << input_path << " contains no numeric rows." << endl;
        return false;
    }

    long rows = temp_rows.size();
    long cols = temp_rows[0].size();
    B.SetDims(rows, cols);
    for (long i = 0; i < rows; ++i)
        for (long j = 0; j < cols; ++j)
            B[i][j] = temp_rows[i][j];

    return true;
}

template <typename T>
bool write_basis(const Mat<T>& M, const string& filename) {
    ofstream out(filename);
    if (!out) {
        cerr << "Failed to open output file: " << filename << endl;
        return false;
    }
    out << M << endl;
    out.close();
    cout << "Matrix written to " << filename << endl;
    return true;
}


template<typename T>
void print_vec(const T& v) {
    cout << "[";
    for (long i = 0; i < v.length(); ++i) {
        cout << v[i];
        if (i < v.length() - 1) cout << " ";
    }
    cout << "]" << endl;
}

template<typename T>
void print_vector_vec(const list<T>& vecs) {
    for (const auto& v : vecs) {
        print_vec(v);
    }
}

template<typename T>
void print_stack_vec(stack<T> vecs) {
    while (!vecs.empty()) {
        print_vec(vecs.top());
        vecs.pop();
    }
}

template<typename T>
bool write_vec_csv(const T& v, const string& filename) {
    ofstream ofs(filename);
    if (!ofs) return false;
    ofs << "[";
    for (long i = 0; i < v.length(); ++i) {
        ofs << v[i];
        if (i < v.length() - 1) ofs << " ";
    }
    ofs << "]" << endl;
    ofs.close();
    return true;
}

template<typename T>
bool write_vector_vec_csv(const list<T>& vecs, const string& filename) {
    ofstream ofs(filename);
    if (!ofs) return false;
    for (const auto& v : vecs) {
        ofs << "[";
        for (long i = 0; i < v.length(); ++i) {
            ofs << v[i];
            if (i < v.length() - 1) ofs << " ";
        }
        ofs << "]" << endl;
    }
    ofs.close();
    return true;
}

template<typename T>
bool write_stack_vec_csv(const stack<T>& vecs, const string& filename) {
    ofstream ofs(filename);
    if (!ofs) return false;
    stack<T> temp_stack = vecs;
    while (!temp_stack.empty()) {
        const T& v = temp_stack.top();
        ofs << "[";
        for (long i = 0; i < v.length(); ++i) {
            ofs << v[i];
            if (i < v.length() - 1) ofs << " ";
        }
        ofs << "]" << endl;
        temp_stack.pop();
    }
    ofs.close();
    return true;
}

#endif // TOOL_H
