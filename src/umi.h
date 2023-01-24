#ifndef _UMI_H_
#define _UMI_H_
#include <string>
#include <vector>
#include "interval.h"
#include "mdf.h"

void add_UMIs(std::vector<pcr_copy> &copies, ostream &umifile, const std::string &format, const std::string &format_back ="");

#endif
