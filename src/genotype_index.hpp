#ifndef GENOTYPE_INDEX_HPP
#define GENOTYPE_INDEX_HPP

#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

// A genotype index from R: a non-negative decimal integer, possibly left-padded by format()
inline uint64_t parse_genotype_index(const std::string & t_index)
{
  const char* start = t_index.c_str();
  while(std::isspace(static_cast<unsigned char>(*start))) start++;
  char* end;
  errno = 0;
  unsigned long long index = std::strtoull(start, &end, 10);
  if(!std::isdigit(static_cast<unsigned char>(*start)) || *end != '\0' || errno == ERANGE){
    throw std::invalid_argument("Invalid genotype index '" + t_index + "'");
  }
  return index;
}

// Indices of marker t_i and of the marker read before it (0 for the first marker and for VCF)
inline void parse_genotype_indices(const std::string & t_genoType,
                                   const std::vector<std::string> & t_genoIndex_prev,
                                   const std::vector<std::string> & t_genoIndex,
                                   size_t t_i,
                                   uint64_t & t_gIndex_prev,
                                   uint64_t & t_gIndex)
{
  t_gIndex = parse_genotype_index(t_genoIndex.at(t_i));
  t_gIndex_prev = 0;
  if(t_i > 0 && t_genoType == "bgen"){
    t_gIndex_prev = parse_genotype_index(t_genoIndex_prev.at(t_i - 1));
  }else if(t_i > 0 && t_genoType == "plink"){
    t_gIndex_prev = parse_genotype_index(t_genoIndex.at(t_i - 1));
  }
}

#endif
