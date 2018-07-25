#ifndef CLSKMERALGORITHM_H
#define CLSKMERALGORITHM_H

#include "string"
#include "stdint.h" // uint64_t
using namespace std;

unsigned int ConvertKmerToNum32(string strKmer);  // convert kmer to number (int) --> the maximum kmer length is 16
uint64_t ConvertKmerToNum64(string strKmer);    // Use int 64 do the things like this

string ConvertNum64ToKmer(uint64_t i64Value, int iKmerLen); // Convert uint64 to string
string ConvertNum32ToKmer(unsigned int i32Value, int iKmerLen); // Convert unsigned int to string

#endif // CLSKMERALGORITHM_H
