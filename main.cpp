//
//  main.cpp
//  Palgirism_Checker
//
//  Created by Omar shaalan on 5/12/21.
//  Copyright © 2021 Omar shaalan. All rights reserved.
//

#include <iostream>
using namespace std;
#include <vector>
#include <fstream>
# define NO_OF_CHARS 256
#include <ctime>

string filetostring(string fileName);
//Brute Force and Hamming Distance
void searchBFHD(string pat, string txt);

//Rabin Karp
void searchRabinKarp(string pat, string text);

//KMP
void computeLPSArray(string pat, int M, int* lps);
void searchKMP(string pat, string txt);

//Boyer Moore Bad Character
int max(int n1, int n2);
void badCharArray(string str, int size, int badchar[NO_OF_CHARS]);
void searchBoyerMooreBadChar(string pattern, string text);

//Boyer Moore Good Suffix
void preprocessStrongSuffixCase1(int* shift, int* bpos, string pat, int m);
void preprocessStrongSuffixCase2(int* shift, int* bpos, string, int m);
void searchBoyerMooreGoodSuffix(string pat, string text);

int main()
{

    //reading the text and pattern files
    string text = filetostring("text20.txt");
    //cout << text<<endl;
    string pattern = filetostring("pattern20.txt");
    //cout << pattern << endl;
    string temp;
    string temppattern = pattern;
    vector <string> patterns;
    //splitting the pattern file to sentences to treat each one as a potential pattern
    for (int i = 0; i < pattern.size() - 1; i++)
        if (pattern[i] == '.')
        {
            temp = pattern.substr(0, i);
            patterns.push_back(temp);
            pattern = pattern.substr(i + 2, pattern.size() - i - 2);
            i = 0;
        }
    pattern = pattern.substr(0, pattern.size() - 1);
    patterns.push_back(pattern);
    pattern = temppattern;
    //calling the four algorithms
    //algorithm 1

    cout << "Brute Force and Hamming Distance\n";
    int start_time1 = clock();
    for (int i = 0; i < patterns.size(); i++)
    {
        searchBFHD((patterns[i]), text);
        cout << ": " << (patterns[i]) << endl;
    }
    int end_time1 = clock();
    cout << endl;
    double execution_time = (double)(end_time1 - start_time1) / (CLOCKS_PER_SEC);
    //execution_time = difftime(end_time1, start_time1);
    cout << execution_time << endl;
    //algorithm 2
    cout << "Rabin Karp\n";
    int start_time2 = clock();
    for (int i = 0; i < patterns.size(); i++)
    {
        searchRabinKarp((patterns[i]), text);
        cout << ": " << (patterns[i]) << endl;
    }
    int end_time2 = clock();

    cout << endl;
    execution_time = (double)(end_time2 - start_time2) / (CLOCKS_PER_SEC);
    //execution_time = difftime(end_time2, start_time2);
    cout << execution_time << endl;
    //algorithm 3
    cout << "KMP\n";
    int start_time3 = clock();
    for (int i = 0; i < patterns.size(); i++)
    {
        searchKMP((patterns[i]), text);
        cout << ": " << (patterns[i]) << endl;
    }
    int end_time3 = clock();
    execution_time = (double)(end_time3 - start_time3) / (CLOCKS_PER_SEC);
    //execution_time = difftime(end_time3, start_time3);
    cout << execution_time << endl;
    cout << endl;
    //algorithm 4
    cout << "Boyer Moore: Bad Character approach\n";
    int start_time4 = clock();
    for (int i = 0; i < patterns.size(); i++)
    {
        searchBoyerMooreBadChar((patterns[i]), text);
        cout << ": " << (patterns[i]) << endl;
    }
    int end_time4 = clock();
    cout << endl;
    execution_time = (double)(end_time4 - start_time4) / (CLOCKS_PER_SEC);
    //execution_time = difftime(end_time4, start_time4);
    cout << execution_time << endl;
    cout << "Boyer Moore: Good Suffix approach\n";
    int start_time5 = clock();
    for (int i = 0; i < patterns.size(); i++)
    {
        searchBoyerMooreGoodSuffix((patterns[i]), text);
        cout << ": " << (patterns[i]) << endl;
    }
    int end_time5 = clock();
    cout << endl;
    execution_time = (double)(end_time5 - start_time5) / (CLOCKS_PER_SEC);
    //execution_time = difftime(end_time5, start_time5);
    cout << execution_time << endl;
    // call an algorithm

    //double execution_time = (double)(end_time - start_time) / (CLOCKS_PER_SEC);
    //cout << execution_time << endl;
    return 0;
}
string filetostring(string fileName) {
    ifstream file(fileName, ios::binary);
    string fileStr;

    istreambuf_iterator<char> inputIt(file), emptyInputIt;
    back_insert_iterator<string> stringInsert(fileStr);

    copy(inputIt, emptyInputIt, stringInsert);

    return fileStr;
}

void searchBFHD(string pat, string txt)
{
    int M = pat.size();
    int N = txt.size();
    int maxMissMatch = 0;
    for (int i = 0; i <= N - M; i++)
    {
        int nmm = 0;
        int j;

        for (j = 0; j < M; j++)
            if (txt[i + j] != pat[j])
            {
                nmm += 1;
                if (nmm > maxMissMatch)
                    break;
            }
        if (j == M)
            cout << "Plagirism found at index " << i;
    }
}

void searchRabinKarp(string pat, string text)
{
    int M = pat.size();
    int N = text.size();
    int i, j;
    int p = 0;
    int t = 0;
    int h = 1;
    int Prime_Number = 101;
    for (i = 0; i < M - 1; i++)
        h = (h * NO_OF_CHARS) % Prime_Number;

    for (i = 0; i < M; i++)
    {
        p = (NO_OF_CHARS * p + pat[i]) % Prime_Number;
        t = (NO_OF_CHARS * t + text[i]) % Prime_Number;
    }

    for (i = 0; i <= N - M; i++)
    {


        if (p == t)
        {
            for (j = 0; j < M; j++)
            {
                if (text[i + j] != pat[j])
                    break;
            }

            if (j == M)
                cout << "Plagirism found at index " << i;
        }


        if (i < N - M) //
        {

            t = (NO_OF_CHARS * (t - text[i] * h) + text[i + M]) % Prime_Number;


            if (t < 0)
                t = (t + Prime_Number);
        }
    }
}

void computeLPSArray(string pat, int M, int* lps)
{
    int len = 0;

    lps[0] = 0;

    int i = 1;
    while (i < M) {
        if (pat[i] == pat[len]) {
            len++;
            lps[i] = len;
            i++;
        }
        else
        {

            if (len != 0) {
                len = lps[len - 1];

            }
            else
            {
                lps[i] = 0;
                i++;
            }
        }
    }
}
void searchKMP(string pat, string txt)
{
    int M = pat.size();
    int N = txt.size();
    int* lps = new int[M];
    computeLPSArray(pat, M, lps);

    int i = 0;
    int j = 0;
    while (i < N) {
        if (pat[j] == txt[i]) {
            j++;
            i++;
        }

        if (j == M) {
            cout << "Plagirism found at index " << i - j;
            j = lps[j - 1];
        }


        else if (i < N && pat[j] != txt[i]) {
            if (j != 0)
                j = lps[j - 1];
            else
                i = i + 1;
        }
    }
}

int max(int n1, int n2) {
    return (n1 > n2) ? n1 : n2;
}
void badCharArray(string pattern, int size, int badchar[NO_OF_CHARS])
{
    int i;

    for (i = 0; i < NO_OF_CHARS; i++)
        badchar[i] = -1;

    for (i = 0; i < size; i++)
        badchar[(int)pattern[i]] = i;
}
void searchBoyerMooreBadChar(string pattern, string text)
{
    int m = pattern.size();
    int n = text.size();

    int badchar[NO_OF_CHARS];

    badCharArray(pattern, m, badchar);

    int nOfShifts = 0;
    while (nOfShifts <= (n - m))
    {
        int j = m - 1;

        while (j >= 0 && pattern[j] == text[nOfShifts + j])
            j--;

        if (j == -1)
        {
            cout << "Plagirism found at index " << nOfShifts;
            if (nOfShifts + m < n)
                nOfShifts = nOfShifts + m - badchar[text[nOfShifts + m]];
            else
                nOfShifts++;
            // nOfShifts+= (nOfShifts+ m < n) ? m - badchar[text[nOfShifts+ m]] : 1;

        }

        else
            nOfShifts += max(1, j - badchar[text[nOfShifts + j]]);
    }
}
void preprocessStrongSuffixCase1(int* shift, int* bpos, string pat, int m)
{
    int i = m, j = m + 1;
    bpos[i] = j;

    while (i > 0)
    {

        while (j <= m && pat[i - 1] != pat[j - 1])
        {
            if (shift[j] == 0)
                shift[j] = j - i;

            j = bpos[j];
        }

        i--; j--;
        bpos[i] = j;
    }
}
void preprocessStrongSuffixCase2(int* shift, int* bpos, string, int m)
{
    int i, j;
    j = bpos[0];
    for (i = 0; i <= m; i++)
    {

        if (shift[i] == 0)
            shift[i] = j;

        if (i == j)
            j = bpos[j];
    }
}
void searchBoyerMooreGoodSuffix(string pat, string text)
{
    int s = 0, j;
    int m = pat.size();
    int n = text.size();
    int* bpos = new int[m + 1];
    int* shift = new int[m + 1];

    for (int i = 0; i < m + 1; i++) shift[i] = 0;

    preprocessStrongSuffixCase1(shift, bpos, pat, m);
    preprocessStrongSuffixCase2(shift, bpos, pat, m);

    while (s <= n - m)
    {

        j = m - 1;

        while (j >= 0 && pat[j] == text[s + j])
            j--;

        if (j < 0)
        {
            cout << "Plagirism found at index " << s;
            s += shift[0];
        }
        else

            s += shift[j + 1];
    }
}

